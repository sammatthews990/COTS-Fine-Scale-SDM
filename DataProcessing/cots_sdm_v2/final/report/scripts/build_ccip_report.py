"""Assemble rendered Quarto content inside the retained CCIP DOCX package.

Only the working output is edited. The source template's styles, page furniture,
opaque parts and original media are retained; instruction/example pages are
replaced with the report. Native Word subsequently refreshes live fields.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import posixpath,math
import re
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from xml.dom import minidom as D

W="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
REL="http://schemas.openxmlformats.org/package/2006/relationships"
CT="http://schemas.openxmlformats.org/package/2006/content-types"
TEMPLATE_SHA="e88e42c58b01be83287e0848dc4ff68b3f5f4f0e3ccc27ba08793424f685a1f5"

def nodes(e,tag): return list(e.getElementsByTagName(tag))
def children(e,tag=None):
    return [n for n in e.childNodes if n.nodeType==n.ELEMENT_NODE and (tag is None or n.tagName==tag)]
def first(e,tag):
    return next((n for n in children(e,tag)),None)
def text(e): return ''.join(n.firstChild.data for n in nodes(e,'w:t') if n.firstChild)
def set_text_node(e,t):
    for c in list(e.childNodes):e.removeChild(c)
    e.appendChild(e.ownerDocument.createTextNode(t))
    e.setAttribute('xml:space','preserve')
def set_slot(e,t):
    ts=nodes(e,'w:t')
    if not ts:raise ValueError('Slot has no text')
    set_text_node(ts[0],t)
    for n in ts[1:]:set_text_node(n,'')
def remove(e):
    if e.parentNode is not None:e.parentNode.removeChild(e)
def el(doc,tag,attrs=None):
    ns=W if tag.startswith('w:') else None
    n=doc.createElementNS(ns,tag)
    for k,v in (attrs or {}).items():n.setAttribute(k,str(v))
    return n
def run(doc,t,style=None):
    r=el(doc,'w:r')
    if style:
        pr=el(doc,'w:rPr');pr.appendChild(el(doc,'w:rStyle',{'w:val':style}));r.appendChild(pr)
    n=el(doc,'w:t');set_text_node(n,t);r.appendChild(n);return r
def field_runs(doc,instruction,result,style=None):
    out=[]
    for kind in ['begin','instruction','separate','result','end']:
        if kind=='result':out.append(run(doc,result,style));continue
        r=el(doc,'w:r')
        if kind=='instruction':
            n=el(doc,'w:instrText',{'xml:space':'preserve'});n.appendChild(doc.createTextNode(' '+instruction+' '))
        else:n=el(doc,'w:fldChar',{'w:fldCharType':kind})
        r.appendChild(n);out.append(r)
    return out
def style_of(p):
    # Quarto can emit two pPr nodes for a caption; Word uses the later style.
    found=[s for pr in children(p,'w:pPr') for s in children(pr,'w:pStyle')]
    return found[-1].getAttribute('w:val') if found else ''
def set_style(p,style):
    props=children(p,'w:pPr');pr=props[0] if props else None
    for extra in props[1:]:
        for child in list(children(extra)):
            for previous in children(pr,child.tagName):remove(previous)
            pr.appendChild(child)
        remove(extra)
    if pr is None:pr=el(p.ownerDocument,'w:pPr');p.insertBefore(pr,p.firstChild)
    s=first(pr,'w:pStyle')
    if s is None:s=el(p.ownerDocument,'w:pStyle');pr.insertBefore(s,pr.firstChild)
    s.setAttribute('w:val',style)
    return pr
def para(doc,style,t=''):
    p=el(doc,'w:p');set_style(p,style)
    if t:p.appendChild(run(doc,t))
    return p
def clone_into(doc,e): return doc.importNode(e,True)
def section_paragraph(doc,section):
    p=para(doc,'BodyText1');first(p,'w:pPr').appendChild(clone_into(doc,section));return p
def page_break(doc):
    p=para(doc,'BodyText1');r=el(doc,'w:r');r.appendChild(el(doc,'w:br',{'w:type':'page'}));p.appendChild(r);return p
def strip_annotations(e):
    for tag in ['w:commentRangeStart','w:commentRangeEnd','w:commentReference']:
        for n in nodes(e,tag):remove(n)

def build(template:Path,content:Path,output:Path,audit:Path,metadata:Path|None=None):
    digest=hashlib.sha256(template.read_bytes()).hexdigest()
    if digest!=TEMPLATE_SHA:raise ValueError('Template changed; inspect and update CCIP_TEMPLATE_GUIDE.md before rebuilding.')
    with ZipFile(template) as z:parts={n:z.read(n) for n in z.namelist()}
    original=dict(parts)
    with ZipFile(content) as z:raw={n:z.read(n) for n in z.namelist()}
    doc=D.parseString(parts['word/document.xml']);body=nodes(doc,'w:body')[0]
    source_nodes=children(body)
    raw_doc=D.parseString(raw['word/document.xml']);raw_body=nodes(raw_doc,'w:body')[0]
    meta=json.loads(metadata.read_text(encoding='utf-8')) if metadata else {}
    raw_paras=nodes(raw_body,'w:p')
    report_title=next(text(p) for p in raw_paras if style_of(p)=='Title')
    report_authors='; '.join(text(p) for p in raw_paras if style_of(p)=='Author') or '[TO CONFIRM: author names and order]'
    report_date=next((text(p) for p in raw_paras if style_of(p)=='Date'),'2026')
    report_year=re.search(r'\b(20\d{2})\b',report_date).group(1)
    for i in range(raw_doc.documentElement.attributes.length):
        a=raw_doc.documentElement.attributes.item(i)
        if a.name.startswith('xmlns:') and not doc.documentElement.hasAttribute(a.name):doc.documentElement.setAttribute(a.name,a.value)

    # Append only missing styles. Existing template definitions win verbatim.
    styles=D.parseString(parts['word/styles.xml']);raw_styles=D.parseString(raw['word/styles.xml'])
    known={x.getAttribute('w:styleId') for x in nodes(styles,'w:style')}
    added_styles=[]
    for st in children(raw_styles.documentElement,'w:style'):
        name=st.getAttribute('w:styleId')
        if name not in known:
            styles.documentElement.appendChild(clone_into(styles,st));known.add(name);added_styles.append(name)
    parts['word/styles.xml']=styles.toxml(encoding='UTF-8')

    rels=D.parseString(parts['word/_rels/document.xml.rels']);raw_rels=D.parseString(raw['word/_rels/document.xml.rels'])
    maxrid=max(int(x.getAttribute('Id')[3:]) for x in children(rels.documentElement) if re.match(r'rId\d+$',x.getAttribute('Id')))
    rid_map={};new_parts=[]
    types=D.parseString(parts['[Content_Types].xml'])
    type_extensions={x.getAttribute('Extension') for x in children(types.documentElement) if x.tagName=='Default'}
    raw_types=D.parseString(raw['[Content_Types].xml'])
    for rel in children(raw_rels.documentElement):
        typ=rel.getAttribute('Type');target=rel.getAttribute('Target')
        if typ.endswith('/image') or typ.endswith('/hyperlink'):
            maxrid+=1;newid=f'rId{maxrid}';rid_map[rel.getAttribute('Id')]=newid
            nr=clone_into(rels,rel);nr.setAttribute('Id',newid)
            if typ.endswith('/image'):
                src=posixpath.normpath(posixpath.join('word',target))
                dest='word/media/report_'+posixpath.basename(target)
                if dest in parts and parts[dest]!=raw[src]:raise ValueError('Media collision '+dest)
                parts[dest]=raw[src];new_parts.append(dest);nr.setAttribute('Target',posixpath.relpath(dest,'word'))
                ext=Path(dest).suffix.lstrip('.')
                if ext not in type_extensions:
                    ct=next(x for x in children(raw_types.documentElement) if x.tagName=='Default' and x.getAttribute('Extension')==ext)
                    types.documentElement.appendChild(clone_into(types,ct));type_extensions.add(ext)
            rels.documentElement.appendChild(nr)

    numbering=D.parseString(parts['word/numbering.xml']);rn=D.parseString(raw['word/numbering.xml'])
    nums={x.getAttribute('w:numId'):x for x in children(numbering.documentElement,'w:num')}
    raw_nums={x.getAttribute('w:numId'):x for x in children(rn.documentElement,'w:num')}
    abstracts={x.getAttribute('w:abstractNumId'):x for x in children(rn.documentElement,'w:abstractNum')}
    numbered_abstract=nodes(nums['4'],'w:abstractNumId')[0].getAttribute('w:val')
    maxnum=max(map(int,nums));ordered_map={}
    def list_type(numid,level):
        a=nodes(raw_nums[numid],'w:abstractNumId')[0].getAttribute('w:val')
        levels=nodes(abstracts[a],'w:lvl')
        lv=next((x for x in levels if x.getAttribute('w:ilvl')==str(level)),levels[0])
        return nodes(lv,'w:numFmt')[0].getAttribute('w:val')
    def number_id(oldid):
        nonlocal maxnum
        if oldid not in ordered_map:
            maxnum+=1;n=el(numbering,'w:num',{'w:numId':maxnum});n.appendChild(el(numbering,'w:abstractNumId',{'w:val':numbered_abstract}))
            ov=el(numbering,'w:lvlOverride',{'w:ilvl':0});ov.appendChild(el(numbering,'w:startOverride',{'w:val':1}));n.appendChild(ov)
            numbering.documentElement.appendChild(n);ordered_map[oldid]=str(maxnum)
        return ordered_map[oldid]

    # Source component patterns, including its explicit portrait/landscape
    # section definitions. Main section numbering restarts at Arabic 1.
    main_section=nodes(source_nodes[136],'w:sectPr')[0]
    portrait=nodes(source_nodes[183],'w:sectPr')[0]
    landscape=nodes(source_nodes[180],'w:sectPr')[0]
    table_pattern=source_nodes[163]
    tbl_pr=first(table_pattern,'w:tblPr')
    body_ppr=first(source_nodes[70],'w:pPr')
    body_ppr=body_ppr.cloneNode(True)
    header_pprs={
        'Heading1':first(source_nodes[69],'w:pPr'),
        'Heading2':first(source_nodes[139],'w:pPr'),
        'Heading3':next(first(p,'w:pPr') for p in source_nodes if p.tagName=='w:p' and style_of(p)=='Heading3'),
        'ReferencesandAppendix':first(source_nodes[65],'w:pPr')}
    source_bullet=next(p for p in source_nodes if p.tagName=='w:p' and style_of(p)=='ListBullet')
    maxbookmark=max([int(n.getAttribute('w:id')) for n in nodes(doc,'w:bookmarkStart')]+[0])+100
    bookmark_ids={};caption_counts={'Figure':0,'Table':0};caption_bookmarks={}
    def normalize_paragraph(p,in_table=False,is_header=False):
        nonlocal maxbookmark
        oldstyle=style_of(p)
        oldpr=first(p,'w:pPr');oldnum=first(oldpr,'w:numPr') if oldpr else None
        heading=oldstyle in ('Heading1','Heading2','Heading3','Heading4','Heading5','Heading6')
        title=text(p).strip()
        caption=bool(re.match(r'^(Figure|Table)\s+\d+[:.]',title)) and ('Caption' in oldstyle)
        picture=bool(nodes(p,'w:drawing'))
        if caption:newstyle='Caption'
        elif picture:newstyle='Normal'
        elif in_table:newstyle='TableHeading' if is_header else 'TableText'
        elif title=='Project Summary' or re.match(r'^Appendix [A-Z]:',title):newstyle='ReferencesandAppendix'
        elif heading:newstyle=oldstyle if oldstyle in header_pprs else 'Heading3'
        elif oldstyle in ['Bibliography','CSLParagraph']:newstyle=oldstyle
        else:newstyle='BodyText1'
        pp=set_style(p,newstyle)
        # Remove imported presentation overrides so the template actually
        # controls body spacing/type, while retaining links and mathematical text.
        for tag in ['w:spacing','w:ind','w:jc','w:tabs','w:shd','w:rPr','w:sectPr']:
            for x in children(pp,tag):remove(x)
        if heading or newstyle=='ReferencesandAppendix':
            pattern=header_pprs[newstyle]
            for c in children(pattern):
                if c.tagName not in ('w:pStyle','w:sectPr','w:rPr'):pp.appendChild(clone_into(doc,c))
        for rpr in nodes(p,'w:rPr'):
            # Keep emphasis, links, scripts and native equation properties.
            for tag in ['w:rFonts','w:sz','w:szCs','w:color','w:highlight']:
                for x in children(rpr,tag):remove(x)
        if nodes(p,'m:oMath') or nodes(p,'m:oMathPara'):
            pp.appendChild(el(doc,'w:spacing',{'w:line':280,'w:lineRule':'atLeast'}))
        if picture:
            pp.appendChild(el(doc,'w:jc',{'w:val':'center'}))
            pp.appendChild(el(doc,'w:keepNext'))
            for ext in nodes(p,'wp:extent'):
                cx=int(ext.getAttribute('cx'));cy=int(ext.getAttribute('cy'));limit=9064*635
                if cx>limit:
                    ratio=limit/cx;ext.setAttribute('cx',str(limit));ext.setAttribute('cy',str(round(cy*ratio)))
                    for ae in nodes(p,'a:ext'):
                        if ae.hasAttribute('cx') and ae.hasAttribute('cy'):
                            ae.setAttribute('cx',str(round(int(ae.getAttribute('cx'))*ratio)))
                            ae.setAttribute('cy',str(round(int(ae.getAttribute('cy'))*ratio)))
        if oldnum and not heading:
            num=first(oldnum,'w:numId');lv=first(oldnum,'w:ilvl')
            oid=num.getAttribute('w:val');level=int(lv.getAttribute('w:val')) if lv else 0
            bullet=list_type(oid,level)=='bullet'
            set_style(p,('ListBullet2' if level else 'ListBullet') if bullet else 'NumberedList')
            for n in children(pp,'w:numPr'):remove(n)
            np=el(doc,'w:numPr');np.appendChild(el(doc,'w:ilvl',{'w:val':min(level,1) if bullet else 0}))
            np.appendChild(el(doc,'w:numId',{'w:val':'2' if bullet else number_id(oid)}));pp.appendChild(np)
        if caption:
            match=re.match(r'^(Figure|Table)\s+(\d+)[:.]\s*',text(p))
            kind,num=match.group(1),int(match.group(2));caption_counts[kind]+=1
            if num!=caption_counts[kind]:raise ValueError(f'{kind} caption order: {num} != {caption_counts[kind]}')
            remaining=match.end()
            for tn in nodes(p,'w:t'):
                val=tn.firstChild.data if tn.firstChild else ''
                take=min(remaining,len(val));set_text_node(tn,val[take:]);remaining-=take
                if remaining==0:break
            name=f'ccip{kind}{num}';caption_bookmarks[(kind,num)]=name
            maxbookmark+=1;bid=str(maxbookmark)
            prefix=[run(doc,kind+' '),el(doc,'w:bookmarkStart',{'w:id':bid,'w:name':name})]
            prefix+=field_runs(doc,f'SEQ {kind} \\* ARABIC',str(num))
            prefix+=[el(doc,'w:bookmarkEnd',{'w:id':bid}),run(doc,': ')]
            before=next((x for x in p.childNodes if x is not pp),None)
            for x in prefix:p.insertBefore(x,before)
            pp.appendChild(el(doc,'w:keepLines'))
            if kind=='Table':pp.appendChild(el(doc,'w:keepNext'))
        return p

    def flatten(e):
        # Pandoc places captions/figures in a one-cell wrapper table. The
        # template uses ordinary picture/caption paragraphs and real data tables.
        if e.tagName=='w:tbl':
            rows=children(e,'w:tr')
            if len(rows)==1 and len(children(rows[0],'w:tc'))==1:
                cell=children(rows[0],'w:tc')[0]
                if any('Caption' in style_of(p) for p in children(cell,'w:p')):
                    return [x for ch in children(cell) if ch.tagName!='w:tcPr' for x in flatten(ch)]
        return [e]
    content_nodes=[]
    started=False
    for e in children(raw_body):
        if e.tagName=='w:p' and text(e).strip()=='Project Summary':started=True
        if not started or e.tagName=='w:sectPr':continue
        if e.tagName=='w:p' and text(e).strip()=='Appendices':continue
        content_nodes.extend(flatten(clone_into(doc,e)))
    if not content_nodes:raise ValueError('Missing Project Summary')

    # Rebase relationships, bookmark IDs and drawing IDs before insertion.
    drawing_id=1000
    for e in content_nodes:
        strip_annotations(e)
        for n in [e]+list(e.getElementsByTagName('*')):
            if n.nodeType!=n.ELEMENT_NODE:continue
            for i in range(n.attributes.length):
                a=n.attributes.item(i)
                if a.namespaceURI==R and a.value in rid_map:a.value=rid_map[a.value]
            if n.tagName in ('w:bookmarkStart','w:bookmarkEnd'):
                old=n.getAttribute('w:id')
                if old not in bookmark_ids:maxbookmark+=1;bookmark_ids[old]=str(maxbookmark)
                n.setAttribute('w:id',bookmark_ids[old])
            if n.tagName=='wp:docPr':drawing_id+=1;n.setAttribute('id',str(drawing_id))
        if e.tagName=='w:p':normalize_paragraph(e)
        elif e.tagName=='w:tbl':
            old=first(e,'w:tblPr')
            if old:remove(old)
            e.insertBefore(clone_into(doc,tbl_pr),e.firstChild)
            rows=children(e,'w:tr')
            operational_table='All recorded dives avoided' in text(e)
            coefficient_table=('Coefficient' in text(rows[0]) and 'Symbol' in text(rows[0]))
            ranking_comparison_table=(
                ('Candidate sites selected' in text(rows[0]) and
                 'Change in dive avoidance' in text(rows[0])) or
                '>0.04 retained: same dive count' in text(rows[0]) or
                '0.02–0.04 retained' in text(rows[0]))
            preferred_widths=[4200,2432,2432] if operational_table else None
            if preferred_widths:
                layout=first(first(e,'w:tblPr'),'w:tblLayout')
                if layout is None:layout=el(doc,'w:tblLayout');first(e,'w:tblPr').appendChild(layout)
                layout.setAttribute('w:type','fixed')
            cols=max(sum(int(first(first(tc,'w:tcPr'),'w:gridSpan').getAttribute('w:val')) if first(tc,'w:tcPr') and first(first(tc,'w:tcPr'),'w:gridSpan') else 1 for tc in children(row,'w:tc')) for row in rows)
            total=13438 if cols>=6 else 9064
            tw=first(first(e,'w:tblPr'),'w:tblW')
            if tw is None:tw=el(doc,'w:tblW');first(e,'w:tblPr').appendChild(tw)
            tw.setAttribute('w:w',str(total));tw.setAttribute('w:type','dxa')
            grid=first(e,'w:tblGrid')
            if grid:
                widths=[int(x.getAttribute('w:w')) for x in children(grid,'w:gridCol')]
                total=13438 if cols>=6 else 9064
                scale=total/sum(widths)
                for x,w in zip(children(grid,'w:gridCol'),widths):x.setAttribute('w:w',str(round(w*scale)))
                if preferred_widths:
                    for x,w in zip(children(grid,'w:gridCol'),preferred_widths):x.setAttribute('w:w',str(w))
            for j,row in enumerate(rows):
                trpr=first(row,'w:trPr')
                if trpr is None:trpr=el(doc,'w:trPr');row.insertBefore(trpr,row.firstChild)
                for h in children(trpr,'w:trHeight'):remove(h)
                if j==0 and not children(trpr,'w:tblHeader'):trpr.appendChild(el(doc,'w:tblHeader'))
                for cell_index,tc in enumerate(children(row,'w:tc')):
                    tcp=first(tc,'w:tcPr')
                    if tcp:
                        for h in children(tcp,'w:tcW'):remove(h)
                        if preferred_widths:tcp.appendChild(el(doc,'w:tcW',{'w:w':preferred_widths[cell_index],'w:type':'dxa'}))
                    for p in nodes(tc,'w:p'):normalize_paragraph(p,True,j==0)
            if operational_table or coefficient_table or ranking_comparison_table:
                ps=nodes(e,'w:p')
                for p in ps[:-1]:first(p,'w:pPr').appendChild(el(doc,'w:keepNext'))
                first(ps[-1],'w:pPr').appendChild(el(doc,'w:keepNext',{'w:val':'0'}))

    # Turn numerical cross-references into native REF fields. Section-name and
    # bibliography hyperlinks remain intact.
    for e in content_nodes:
        for h in list(nodes(e,'w:hyperlink')):
            anchor=h.getAttribute('w:anchor');display=text(h)
            kind='Figure' if anchor.startswith('fig') else ('Table' if anchor.startswith('tbl') else None)
            match=re.search(r'(\d+)$',display) if kind else None
            if match and (kind,int(match.group(1))) in caption_bookmarks:
                num=int(match.group(1));prefix=display[:match.start()]
                repl=([run(doc,prefix)] if prefix else [])+field_runs(doc,f'REF {caption_bookmarks[(kind,num)]} \\h',str(num),'Hyperlink')
                for n in repl:h.parentNode.insertBefore(n,h)
                remove(h)

    # Keep an inline image and its caption together within the portrait text
    # height. Tall supplementary plots otherwise create empty/image-only pages.
    for i,e in enumerate(content_nodes):
        if e.tagName!='w:p' or not nodes(e,'w:drawing'):continue
        caption=next((p for p in content_nodes[i+1:i+5] if p.tagName=='w:p' and style_of(p)=='Caption'),None)
        cap_lines=math.ceil(len(text(caption))/95)+1 if caption else 2
        limit=(13438-(cap_lines*240+520))*635
        for ext in nodes(e,'wp:extent'):
            cy=int(ext.getAttribute('cy'));cx=int(ext.getAttribute('cx'))
            if cy>limit:
                ratio=limit/cy;ext.setAttribute('cx',str(round(cx*ratio)));ext.setAttribute('cy',str(limit))
                for ae in nodes(e,'a:ext'):
                    if ae.hasAttribute('cx') and ae.hasAttribute('cy'):
                        ae.setAttribute('cx',str(round(int(ae.getAttribute('cx'))*ratio)))
                        ae.setAttribute('cy',str(round(int(ae.getAttribute('cy'))*ratio)))

    # Express a requested section start on the heading itself. A separate
    # break-only paragraph can overflow and create an otherwise empty page
    # when the template heading also has pageBreakBefore.
    for i in range(len(content_nodes)-2,-1,-1):
        p=content_nodes[i]
        following=next((n for n in content_nodes[i+1:]
            if n.tagName not in ('w:bookmarkStart','w:bookmarkEnd')),None)
        if (p.tagName=='w:p' and not text(p).strip()
            and any(b.getAttribute('w:type')=='page' for b in nodes(p,'w:br'))
            and following is not None and following.tagName=='w:p'
            and style_of(following) in ('Heading1','Heading2','ReferencesandAppendix')):
            pp=first(following,'w:pPr')
            if not first(pp,'w:pageBreakBefore'):pp.appendChild(el(doc,'w:pageBreakBefore'))
            content_nodes.pop(i)

    # Keep the short reproducibility appendix together instead of leaving its
    # final access sentence alone on the following page.
    short_appendix=False;group=[]
    for e in content_nodes:
        if e.tagName=='w:p' and re.match(r'^Appendix [A-Z]:',text(e).strip()):
            if short_appendix:
                for p in group[:-1]:first(p,'w:pPr').appendChild(el(doc,'w:keepNext'))
            short_appendix=text(e).startswith('Appendix F:');group=[]
        if short_appendix and e.tagName=='w:p':group.append(e)

    for n in list(body.childNodes):body.removeChild(n)
    # Preserve source cover/publication components. Text replacements do not
    # touch their original run properties, logos, shapes or section break.
    front=[clone_into(doc,n) for n in source_nodes[:19]]
    cover_rows=children(front[0],'w:tr')
    cover_ps=children(children(cover_rows[0],'w:tc')[0],'w:p')
    set_slot(cover_ps[0],report_title)
    for p in cover_ps[1:]:
        for t in nodes(p,'w:t'):set_text_node(t,'')
    set_slot(cover_rows[1],report_authors)
    affiliation_cell=children(cover_rows[2],'w:tc')[0]
    ap=children(affiliation_cell,'w:p')[0].cloneNode(True)
    for p in list(children(affiliation_cell,'w:p')):remove(p)
    for affiliation in meta.get('affiliations',['[TO CONFIRM: author affiliations]']):
        p=ap.cloneNode(True);set_slot(p,affiliation);affiliation_cell.appendChild(p)
    for idx,val in [(4,meta.get('contact_name','[TO CONFIRM: contact name]')),
        (5,meta.get('contact_organisation','[TO CONFIRM: organisation]')),
        (6,meta.get('contact_email','[TO CONFIRM: contact email]')),
        (8,f'{meta.get("citation_authors","[TO CONFIRM: authors]")} ({report_year}) {report_title}. A report to the Australian Government by the COTS Control Innovation Program ('),
        (17,'© Copyright: '+meta.get('copyright_owner','[TO CONFIRM: copyright owner]')+'.')]:set_slot(front[idx],val)
    for r in field_runs(doc,'NUMPAGES','0'):front[8].appendChild(r)
    front[8].appendChild(run(doc,' pp).'))
    for n in front:strip_annotations(n);body.appendChild(n)

    def heading_from(index):
        p=clone_into(doc,source_nodes[index])
        for tag in ['w:drawing','w:pict']:
            for x in nodes(p,tag):remove(x)
        strip_annotations(p);return p
    def toc_field(index,instruction):
        p=para(doc,style_of(source_nodes[index]) or 'TOC1')
        for r in field_runs(doc,instruction,'Update fields in Word'):p.appendChild(r)
        return p
    body.appendChild(heading_from(19))
    body.appendChild(toc_field(21,'TOC \\o "1-3" \\h \\z \\u'))
    # Put the break on each list heading: a separate break paragraph can
    # create an empty page when the expanded contents fills its last page.
    fig_heading=heading_from(49)
    set_style(fig_heading,style_of(fig_heading)).appendChild(el(doc,'w:pageBreakBefore'))
    body.appendChild(fig_heading)
    body.appendChild(toc_field(51,'TOC \\h \\z \\c "Figure"'))
    table_heading=heading_from(55)
    set_style(table_heading,style_of(table_heading)).appendChild(el(doc,'w:pageBreakBefore'))
    body.appendChild(table_heading)
    body.appendChild(toc_field(57,'TOC \\h \\z \\c "Table"'))
    body.appendChild(clone_into(doc,source_nodes[59]))
    body.appendChild(heading_from(60))
    acronyms=[('COTS','Crown-of-thorns starfish'),('CCIP','COTS Control Innovation Program'),('CPUE','Catch per unit effort'),
      ('GBR','Great Barrier Reef'),('GBRF','Great Barrier Reef Foundation'),('GBRMPA / Reef Authority','Great Barrier Reef Marine Park Authority'),
      ('AIMS','Australian Institute of Marine Science'),('SDM','Species distribution model'),('LTMP','Long-Term Monitoring Program'),
      ('SALAD','Scooter-assisted large area diver-based survey'),('RHIS','Reef Health and Impact Surveys'),('RRAP','Reef Restoration and Adaptation Program'),
      ('NESP','National Environmental Science Program'),('KMZ','Compressed Keyhole Markup Language spatial file'),('FAIR','Findable, accessible, interoperable and reusable')]
    at=clone_into(doc,source_nodes[62]);row=children(at,'w:tr')[0].cloneNode(True)
    for r in list(children(at,'w:tr')):remove(r)
    for a,b in acronyms:
        r=row.cloneNode(True);cells=children(r,'w:tc');set_slot(cells[0],a);set_slot(cells[1],b)
        for p in nodes(r,'w:p'):set_style(p,'TableText')
        at.appendChild(r)
    body.appendChild(at);body.appendChild(clone_into(doc,source_nodes[64]))

    in_appendix=False;wide_tables=0;equation_count=0
    for e in content_nodes:
        if e.tagName=='w:p' and nodes(e,'m:oMathPara'):
            # The supplied template places display equations and numbers in
            # borderless tables. Reuse that component with editable OMML,
            # replacing its legacy Equation.3 objects and example number.
            equation_count+=1
            lead=next((p for p in reversed(children(body)) if p.tagName=='w:p'),None)
            if lead is not None and text(lead).rstrip().endswith(':'):
                first(lead,'w:pPr').appendChild(el(doc,'w:keepNext'))
            eq=clone_into(doc,source_nodes[166]);row=children(eq,'w:tr')[0]
            cells=children(row,'w:tc');left,middle,right=cells
            remove(middle)
            for cell in (left,right):
                for child in list(children(cell)):
                    if child.tagName!='w:tcPr':remove(child)
            lp=first(left,'w:tcPr');first(lp,'w:tcW').setAttribute('w:w','6520')
            lp.appendChild(el(doc,'w:gridSpan',{'w:val':'2'}))
            # Let the native equation expand vertically instead of clipping.
            for h in list(nodes(row,'w:trHeight')):remove(h)
            pp=first(e,'w:pPr');pp.appendChild(el(doc,'w:keepNext',{'w:val':'0'}))
            left.appendChild(e)
            number=para(doc,'BodyText1');np=first(number,'w:pPr')
            np.appendChild(el(doc,'w:jc',{'w:val':'right'}))
            number.appendChild(run(doc,'('))
            for r in field_runs(doc,'SEQ Equation \\* ARABIC',str(equation_count)):number.appendChild(r)
            number.appendChild(run(doc,')'));right.appendChild(number)
            body.appendChild(eq)
            continue
        if e.tagName=='w:p' and re.match(r'^Appendix [A-Z]:',text(e).strip()) and not in_appendix:
            # Avoid a redundant manual page break immediately before the
            # template's own section break.
            while children(body) and children(body)[-1].tagName=='w:p' and nodes(children(body)[-1],'w:br') and not text(children(body)[-1]).strip():remove(children(body)[-1])
            body.appendChild(section_paragraph(doc,main_section));in_appendix=True
        wide=e.tagName=='w:tbl' and max(len(children(r,'w:tc')) for r in children(e,'w:tr'))>=6
        if wide and in_appendix:
            # The closing sectPr describes the section before it. Thus portrait
            # closes before a wide table; landscape closes immediately after it.
            # Include the nearby subsection introduction and caption in the
            # wide section, so a caption is not stranded on a portrait page.
            previous=children(body);cut=None;paragraphs=0
            for j in range(len(previous)-1,-1,-1):
                p=previous[j]
                if p.tagName=='w:tbl':break
                if p.tagName=='w:p':
                    paragraphs+=1
                    if style_of(p) in ('Heading1','Heading2','Heading3','ReferencesandAppendix'):
                        cut=j;break
                    if paragraphs>8:break
            if cut is None:
                cut=next((j for j in range(len(previous)-1,max(-1,len(previous)-5),-1)
                    if previous[j].tagName=='w:p' and style_of(previous[j])=='Caption'),len(previous))
            tail=previous[cut:]
            for p in tail:remove(p)
            body.appendChild(section_paragraph(doc,portrait))
            for p in tail:body.appendChild(p)
            body.appendChild(e)
            body.appendChild(section_paragraph(doc,landscape));wide_tables+=1
        else:body.appendChild(e)
    body.appendChild(section_paragraph(doc,portrait))
    back=[clone_into(doc,n) for n in source_nodes[184:]]
    set_slot(back[-2],meta.get('back_contact','[TO CONFIRM: key contact name, title and email]'))
    for n in back:strip_annotations(n);body.appendChild(n)

    # The source comments contain authoring guidance, not review feedback. Clear
    # those four instructions from the filled copy, retaining the package parts.
    strip_annotations(body)
    for name in list(parts):
        if re.match(r'word/comments(?:Extended|Ids)?\.xml$',name):
            cd=D.parseString(parts[name])
            for ch in list(cd.documentElement.childNodes):cd.documentElement.removeChild(ch)
            parts[name]=cd.toxml(encoding='UTF-8')
    for name in list(parts):
        if re.match(r'word/footer\d+\.xml$',name):
            fd=D.parseString(parts[name]);changed=False
            for t in nodes(fd,'w:t'):
                if t.firstChild and '[Project Code]' in t.firstChild.data:
                    set_text_node(t,t.firstChild.data.replace('[Project Code]','CCIP-BR-13'));changed=True
            if changed:parts[name]=fd.toxml(encoding='UTF-8')
    settings=D.parseString(parts['word/settings.xml'])
    uf=nodes(settings,'w:updateFields')
    if not uf:settings.documentElement.appendChild(el(settings,'w:updateFields',{'w:val':'true'}))
    else:uf[0].setAttribute('w:val','true')
    parts['word/settings.xml']=settings.toxml(encoding='UTF-8')
    parts['word/numbering.xml']=numbering.toxml(encoding='UTF-8')
    parts['word/_rels/document.xml.rels']=rels.toxml(encoding='UTF-8')
    parts['[Content_Types].xml']=types.toxml(encoding='UTF-8')
    parts['word/document.xml']=doc.toxml(encoding='UTF-8')
    output.parent.mkdir(parents=True,exist_ok=True)
    with ZipFile(output,'w',ZIP_DEFLATED) as z:
        for n,data in parts.items():z.writestr(n,data)
    changed=[n for n in original if parts[n]!=original[n]]
    preserved=[n for n in original if parts[n]==original[n]]
    # Opaque source parts and source media are never rewritten during assembly.
    for n in original:
        if n.startswith(('customXml/','word/media/')) or n in ('word/theme/theme1.xml','word/fontTable.xml'):
            assert original[n]==parts[n],n
    if hashlib.sha256(template.read_bytes()).hexdigest()!=TEMPLATE_SHA:raise ValueError('Template was modified')
    report={'template_sha256':TEMPLATE_SHA,'template_parts_preserved':preserved,'parts_changed':changed,
      'parts_added':new_parts,'styles_added':added_styles,'captions':caption_counts,'wide_appendix_tables':wide_tables,'editable_display_equations':equation_count,
      'source_comments_removed_from_body':True,'output':str(output)}
    audit.parent.mkdir(parents=True,exist_ok=True);audit.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['captions','wide_appendix_tables','styles_added','parts_changed']},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--template',type=Path,required=True);ap.add_argument('--content',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--metadata',type=Path)
    a=ap.parse_args();build(a.template,a.content,a.output,a.audit,a.metadata)
