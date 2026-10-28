"""Check template identity, critical styles, fields, scope and stated limits."""
import argparse,hashlib,json,re,posixpath,subprocess
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as E

N={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
W='{'+N['w']+'}'
SHA='e88e42c58b01be83287e0848dc4ff68b3f5f4f0e3ccc27ba08793424f685a1f5'
def load(z,n):return E.fromstring(z.read(n))
def ptext(p):return ''.join(x.text or '' for x in p.findall('.//w:t',N))
def main(a):
 assert hashlib.sha256(a.template.read_bytes()).hexdigest()==SHA,'Source template changed'
 with ZipFile(a.template) as t,ZipFile(a.docx) as z:
  doc=load(z,'word/document.xml');styles=load(z,'word/styles.xml');oldstyles=load(t,'word/styles.xml')
  ss={x.get(W+'styleId'):x for x in styles.findall('w:style',N)}
  os={x.get(W+'styleId'):x for x in oldstyles.findall('w:style',N)}
  critical=['Normal','BodyText','BodyText1','Heading1','Heading2','Heading3','Caption','TableText','TableHeading','Contents','ReferencesandAppendix']
  for sid in critical:
   assert sid in ss,sid
   for path in ['w:rPr/w:rFonts','w:rPr/w:sz','w:rPr/w:color','w:rPr/w:b','w:rPr/w:caps','w:pPr/w:spacing']:
    expected=os[sid].find(path,N)
    if expected is not None:
     actual=ss[sid].find(path,N)
     assert actual is not None,(sid,path,'missing')
     for k,v in expected.attrib.items():assert actual.get(k)==v,(sid,path,k,v,actual.get(k))
  alltext='\n'.join(ptext(p) for p in doc.findall('.//w:p',N))
  for phrase in ['Where and when to search for COTS on a reef','Engagement with Rights-holders and Management Stakeholders','Data Accessibility','5 November 2026','January–February 2027']:
   assert phrase.casefold() in alltext.casefold(),phrase
  for phrase in ['Length guidance:','This is the Heading 2 style','Use this section as a quick reference','Report title: report subtitle as needed']:
   assert phrase not in alltext,phrase
  assert not doc.findall('.//w:commentReference',N),'Template guidance comment anchors remain'
  instr=' '.join(x.text or '' for x in doc.findall('.//w:instrText',N))
  for code in ['TOC','SEQ Figure','SEQ Table','REF ccip']:
   assert code in instr,('Missing live field',code)
  template_doc=load(t,'word/document.xml');template_body=list(template_doc.find('w:body',N))
  rels=load(t,'word/_rels/document.xml.rels');targets={r.get('Id'):r.get('Target') for r in rels}
  drawing_ns={'a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
  image_ids={n.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed')
   for p in template_body[:19]+template_body[184:] for n in p.findall('.//a:blip',drawing_ns)}
  saved_media={hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/')}
  # Word may rename media and purge unused guide screenshots during field
  # refresh. The cover/back-page branding itself must remain byte-identical.
  for rid in image_ids:
   name=posixpath.normpath(posixpath.join('word',targets[rid]))
   assert hashlib.sha256(t.read(name)).hexdigest() in saved_media,('Branding image changed',name)
  assert all(n in z.namelist() for n in t.namelist() if n.startswith('customXml/')),'Opaque source part lost'
  for name in ['Funding Acknowledgement','Traditional Owner Acknowledgement','Disclaimer']:
   assert name.casefold() in alltext.casefold(),name
  assert '[Project Code]' not in '\n'.join(z.read(n).decode('utf-8') for n in z.namelist() if re.match(r'word/footer\d+\.xml$',n))
  pdf_text=subprocess.check_output(['pdftotext','-layout',str(a.docx.with_suffix('.pdf')),'-']).decode('utf-8')
  for phrase in ['Error! Reference source not found','Error! Bookmark not defined','No table of figures entries found','Update fields in Word']:
   assert phrase.casefold() not in pdf_text.casefold(),phrase
  for number,page in enumerate(pdf_text.split('\f')[:-1],1):
   content=re.sub(r'Bridging Project CCIP-BR-13','',page)
   content=re.sub(r'Page\s+(?:\d+|[ivxlcdm]+)','',content,flags=re.I)
   assert content.strip(),f'Empty page {number}; inspect section/page breaks'
 source=a.source.read_text(encoding='utf-8')
 display_equations=len(re.findall(r'(?m)^\$\$\s*$',source))//2
 math_ns={'m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
 assert len(doc.findall('.//m:oMathPara',math_ns))==display_equations,'Editable display equation missing'
 if display_equations:assert 'SEQ Equation' in instr,'Live equation numbering missing'
 # Remove code before counting words, and treat each inline result as one value.
 def body(heading,next_heading):
  return source[source.index(heading):source.index(next_heading)]
 def words(s):
  s=re.sub(r'```.*?```','',s,flags=re.S);s=re.sub(r'`r [^`]+`','RESULT',s)
  s=re.sub(r'<!--.*?-->','',s,flags=re.S);return len(re.findall(r'\b[\w–-]+\b',s))
 summary=words(body('# Project Summary','# Context & Objectives'))
 conclusions=words(body('# Conclusions & Next Steps','# Acknowledgements'))
 assert summary<=360,('Summary word budget',summary)
 assert conclusions<=300,('Conclusions exceed template limit',conclusions)
 results=body('# Results','# Discussion & Outputs')
 displays=len(re.findall(r'(?m)^!\[',results))+results.count('#| label: tbl-')
 assert 2<=displays<=4,('Results display count',displays)
 methods=body('# Methods','# Results');mfig=len(re.findall(r'(?m)^!\[',methods))
 assert 1<=mfig<=3,('Methods figure count',mfig)
 pages=json.loads(a.page_map.read_text(encoding='utf-8-sig'))
 page_citation=re.search(r'Program\s*\((\d+)\s+pp\)',alltext)
 assert page_citation and int(page_citation.group(1))==pages['pages'],'Citation NUMPAGES cache is stale'
 headings=pages['headings']
 def find(title):return next(x for x in headings if x['text'].strip().casefold()==title.casefold())
 ps=find('Project Summary');context=find('Context & Objectives')
 assert ps['page']==1,('Summary must start at Arabic 1',ps)
 assert context['physical_page']-ps['physical_page']==1,('Summary exceeds one page',ps,context)
 context_pages=find('Methods')['physical_page']-context['physical_page']
 discussion_pages=find('Outcomes & Impact Pathway')['physical_page']-find('Discussion & Outputs')['physical_page']
 outcome_pages=find('Conclusions & Next Steps')['physical_page']-find('Outcomes & Impact Pathway')['physical_page']
 warnings=[]
 for section,n,lo,hi in [('Context',context_pages,2,3),('Discussion',discussion_pages,3,4),('Outcomes',outcome_pages,1,2)]:
  if not lo<=n<=hi:warnings.append(f'{section}: {n} pages; guidance {lo}–{hi}')
 result={'template_sha256':SHA,'pages':pages['pages'],'summary_words':summary,'conclusion_words':conclusions,
  'main_results_displays':displays,'methods_figures':mfig,'editable_display_equations':display_equations,'context_pages':context_pages,
  'discussion_pages':discussion_pages,'outcome_pages':outcome_pages,'warnings':warnings,'critical_styles_preserved':critical}
 a.out.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for arg in ['template','docx','source','page-map','out']:p.add_argument('--'+arg,type=Path,required=True)
 main(p.parse_args())
