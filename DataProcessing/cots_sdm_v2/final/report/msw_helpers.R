# Helpers for reporting numbers straight from results tables, so no value in the
# text is typed by hand. Source this in the document's setup chunk.
#
#   res <- tab("forward_validation_summary.csv")
#   auc <- pick(res, score = "habitat", subset = "all", col = "auc")
#   `r num(auc)`                    -> 0.57
#   `r pct(0.384)`                  -> 38%
#   `r ci(0.68, 0.64, 0.72)`        -> 0.68 (95% CI 0.64–0.72)
#   `r int(18400)`                  -> 18,400
#
# pick() stops if the filter does not match exactly one row, so a renamed column
# or a changed table breaks the render instead of silently printing a wrong
# number.

msw_tables_dir <- Sys.getenv("MSW_TABLES", "outputs/tables")

tab <- function(file, dir = msw_tables_dir) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Results table not found: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

pick <- function(df, ..., col) {
  filters <- list(...)
  keep <- rep(TRUE, nrow(df))
  for (nm in names(filters)) {
    if (!nm %in% names(df)) stop("pick(): no column '", nm, "'", call. = FALSE)
    v <- filters[[nm]]
    keep <- keep & if (is.numeric(v)) abs(as.numeric(df[[nm]]) - v) < 1e-9 else df[[nm]] %in% v
  }
  if (!col %in% names(df)) stop("pick(): no column '", col, "'", call. = FALSE)
  if (sum(keep) != 1) stop("pick(): ", sum(keep), " rows match ",
                           paste(names(filters), filters, sep = "=", collapse = ", "), call. = FALSE)
  df[[col]][keep]
}

num <- function(x, digits = 2) formatC(as.numeric(x), format = "f", digits = digits, big.mark = ",")
int <- function(x) formatC(round(as.numeric(x)), format = "d", big.mark = ",")
pct <- function(x, digits = 0, scale = 100) paste0(formatC(scale * as.numeric(x), format = "f", digits = digits), "%")
pts <- function(x, digits = 1) paste(formatC(100 * as.numeric(x), format = "f", digits = digits), "percentage points")
ci <- function(est, lo, hi, digits = 2, f = num, level = "95% CI") {
  sprintf("%s (%s %s–%s)", f(est, digits), level, f(lo, digits), f(hi, digits))
}
rng <- function(lo, hi, digits = 2, f = num) sprintf("%s–%s", f(lo, digits), f(hi, digits))

# --- Title block -------------------------------------------------------------
# Authors, numbered affiliations and correspondence, built from the YAML keys
# msw-authors / msw-affiliations so the same block appears in PDF and DOCX
# (Quarto's own title block drops affiliations in both). Use in an R chunk with
# `#| output: asis` straight after the setup chunk.
#
# msw-authors:
#   - name: Jane Q. Citizen
#     affiliations: [aims, jcu]
#     corresponding: true        # optional
#     email: j.citizen@aims.gov.au
#     equal: true                # optional, shared first authorship (†)
# msw-affiliations:
#   - id: aims
#     department: …              # optional
#     name: Australian Institute of Marine Science
#     city: Townsville
#     region: QLD                # optional
#     country: Australia

msw_title_block <- function(path = knitr::current_input(dir = TRUE)) {
  m <- rmarkdown::yaml_front_matter(path)
  au <- m[["msw-authors"]]
  af <- m[["msw-affiliations"]]
  if (!length(au)) return(invisible(""))
  ids <- vapply(af, function(a) as.character(a$id), "")
  used <- unique(unlist(lapply(au, function(a) a$affiliations)))
  unknown <- setdiff(used, ids)
  if (length(unknown)) stop("msw-authors cite unknown affiliation ids: ", paste(unknown, collapse = ", "), call. = FALSE)
  n <- stats::setNames(seq_along(used), used)
  who <- vapply(au, function(a) {
    marks <- c(n[unlist(a$affiliations)], if (isTRUE(a$corresponding)) "\\*", if (isTRUE(a$equal)) "†")
    paste0(a$name, if (length(marks)) paste0("^", paste(marks, collapse = ","), "^"))
  }, "")
  aff <- vapply(used, function(id) {
    a <- af[[match(id, ids)]]
    parts <- unlist(a[intersect(c("department", "name", "city", "region", "country"), names(a))])
    paste0("^", n[[id]], "^ ", paste(parts, collapse = ", "))
  }, "")
  corr <- Filter(function(a) isTRUE(a$corresponding), au)
  eq <- any(vapply(au, function(a) isTRUE(a$equal), logical(1)))
  out <- c(paste(who, collapse = ", "), "", paste(aff, collapse = "\\\n"), "",
           if (length(corr)) paste0("^\\*^ Correspondence: ", paste(vapply(corr, function(a) paste0(a$name, ", ", a$email), ""), collapse = "; ")),
           if (eq) "\n^†^ These authors contributed equally to this work and share first authorship.")
  paste(out, collapse = "\n")
}

# --- Word counts -------------------------------------------------------------
# Counts words in the .qmd source the way journals count "main text": body
# paragraphs only. Excluded: YAML, code chunks, HTML comments, headings, tables,
# figure lines, captions, div fences, citation keys, and everything from the
# first end-matter heading (statements, references, legends) onward.
# Inline `r ...` expressions count as one word. Expect +/- 2 % against Word.

msw_end_matter <- c("data availability", "code availability", "ethics statement", "author contributions",
                    "funding", "acknowledg", "conflict of interest", "competing interests",
                    "generative ai", "publisher", "supplementary", "references", "figure legends",
                    "figure captions", "extended data", "additional information", "declaration")

msw_sections <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  # Drop YAML front matter.
  if (length(x) && x[1] == "---") {
    end <- which(x == "---")[2]
    x <- x[-seq_len(end)]
  }
  txt <- paste(x, collapse = "\n")
  txt <- gsub("(?s)<!--.*?-->", "", txt, perl = TRUE)
  txt <- gsub("(?s)\n```[^\n]*\n.*?\n```", "\n", txt, perl = TRUE)
  x <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  heading <- grepl("^#{1,6} ", x)
  section <- cumsum(heading)
  titles <- c("(before first heading)", sub("\\s*\\{.*\\}\\s*$", "", sub("^#{1,6} ", "", x[heading])))
  levels <- c(0L, nchar(sub("^(#+).*", "\\1", x[heading])))
  # Front-matter lines (keywords, counts, citation) are never body text.
  front_line <- grepl("^\\*\\*(Keywords|Word count|Main text|Article type|Suggested citation|Report No)", x) |
    grepl("^\\s*\\{\\{<.*>\\}\\}\\s*$", x)
  x[front_line] <- ""
  body <- x
  body[heading | grepl("^\\s*\\|", body) | grepl("^\\s*!\\[", body) | grepl("^\\s*:::", body) |
         grepl("^\\s*: ", body) | grepl("^\\s*#\\|", body)] <- ""
  # raw: everything a reader sees in the section, including table cells (used
  # for per-section limits such as the CCIP form's boxes).
  raw <- x
  raw[heading | grepl("^\\s*:::", raw) | grepl("^\\s*#\\|", raw) | grepl("^\\s*\\|[-: |]+\\|\\s*$", raw)] <- ""
  raw <- gsub("|", " ", raw, fixed = TRUE)
  data.frame(section = section, title = titles[section + 1], level = levels[section + 1], text = body,
             raw = raw, stringsAsFactors = FALSE)
}

msw_count_words <- function(text) {
  text <- gsub("`r [^`]*`", " X ", text)
  text <- gsub("\\[[^]]*@[^]]*\\]", " ", text)          # [@key; @key2]
  text <- gsub("(^|\\s)-?@[A-Za-z0-9_:.-]+", " ", text)  # @key and -@key
  text <- gsub("\\{[#.][^}]*\\}", " ", text)             # {#fig-x} {.class}
  text <- gsub("[*_`#>]", " ", text)
  w <- unlist(strsplit(text, "\\s+"))
  sum(grepl("[[:alnum:]]", w))
}

# Main-text words (Frontiers rule) and, for short-form articles, Methods words.
msw_word_count <- function(path = knitr::current_input(dir = TRUE)) {
  s <- msw_sections(path)
  is_end <- vapply(tolower(s$title), function(t) any(startsWith(t, msw_end_matter)), logical(1))
  first_end <- if (any(is_end & s$level == 1)) min(s$section[is_end & s$level == 1]) else Inf
  is_front <- s$level == 1 & tolower(s$title) %in% c("abstract", "summary", "summary paragraph")
  body <- s[s$section >= 1 & s$section < first_end & !is_front, ]
  methods_top <- unique(body$section[body$level == 1 & grepl("^(online )?methods", tolower(body$title))])
  in_methods <- rep(FALSE, nrow(body))
  for (m in methods_top) {
    nxt <- body$section[body$level == 1 & body$section > m]
    in_methods <- in_methods | (body$section >= m & body$section < if (length(nxt)) min(nxt) else Inf)
  }
  abstract <- s$text[s$section %in% s$section[is_front]]
  list(main = msw_count_words(body$text[!in_methods]), methods = msw_count_words(body$text[in_methods]),
       abstract = msw_count_words(abstract),
       by_section = tapply(seq_len(nrow(s)), s$section, function(i) msw_count_words(s$text[i])),
       titles = tapply(s$title, s$section, `[`, 1))
}
