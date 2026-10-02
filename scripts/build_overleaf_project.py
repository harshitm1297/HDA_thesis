"""Convert the canonical Markdown thesis into structured Overleaf chapter files."""

from __future__ import annotations

import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "THESIS_REPORT.md"
TARGET = PROJECT / "overleaf_thesis"


CITATIONS = {
    "(Webb-Robertson et al., 2015; Lazar et al., 2016; O'Brien et al., 2018; Li et al., 2023)": r"\parencite{webb2015,lazar2016,obrien2018,li2023}",
    "(Webb-Robertson et al., 2015; Lazar et al., 2016; O'Brien et al., 2018; Chion et al., 2022)": r"\parencite{webb2015,lazar2016,obrien2018,chion2022}",
    "(Smyth, 2004; Ritchie et al., 2015)": r"\parencite{smyth2004,ritchie2015}",
    "(Subramanian et al., 2005)": r"\parencite{subramanian2005}",
    "(Khatri et al., 2012)": r"\parencite{khatri2012}",
    "(Monti et al., 2003; Şenbabaoğlu et al., 2014)": r"\parencite{monti2003,senbabaoglu2014}",
    "(Varma and Simon, 2006; Lewis et al., 2023)": r"\parencite{varma2006,lewis2023}",
    "(Davis et al., 2023; Kapoor and Narayanan, 2023)": r"\parencite{davis2023,kapoor2023}",
    "(Zou and Hastie, 2005)": r"\parencite{zou2005}",
    "(Meinshausen and Bühlmann, 2010)": r"\parencite{meinshausen2010}",
    "(Steegen et al., 2016)": r"\parencite{steegen2016}",
    "(Wilkinson et al., 2016)": r"\parencite{wilkinson2016}",
    "(Taylor et al., 2007)": r"\parencite{taylor2007}",
    "(Benjamini and Hochberg, 1995)": r"\parencite{benjamini1995}",
    "(Smyth, 2004; Ritchie et al., 2015)": r"\parencite{smyth2004,ritchie2015}",
    "(Moore et al., 2011)": r"\parencite{moore2011}",
    "(Collins et al., 2015; Wolff et al., 2019)": r"\parencite{collins2015,wolff2019}",
    "(Li and Smyth, 2023)": r"\parencite{li2023}",
    "Li and Smyth (2023)": r"\textcite{li2023}",
    "Regulation (EU) 2017/746": r"\textcite{eu2017ivdr}",
    "Khatri et al. (2012)": r"\textcite{khatri2012}",
    "Davis et al. (2023)": r"\textcite{davis2023}",
    "Kapoor and Narayanan (2023)": r"\textcite{kapoor2023}",
    "(s_{j,*}^2)": r"\(s_{j,*}^{2}\)",
    "(Chion et al., 2022)": r"\parencite{chion2022}",
    "(Varma and Simon, 2006)": r"\parencite{varma2006}",
}


def escape_plain(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def inline(text: str) -> str:
    tokens: list[str] = []

    def protect(value: str) -> str:
        token = f"@@TOKEN{len(tokens):04d}@@"
        tokens.append(value)
        return token

    for source, latex in sorted(CITATIONS.items(), key=lambda item: -len(item[0])):
        text = text.replace(source, protect(latex))

    text = re.sub(
        r"\[([^\]]+)\]\((https?://[^)]+)\)",
        lambda match: protect(r"\href{" + match.group(2) + "}{" + escape_plain(match.group(1)) + "}"),
        text,
    )
    def render_code(match: re.Match[str]) -> str:
        value = match.group(1)
        # URL's \path form permits line breaks at separators in long file paths
        # and underscore-heavy machine identifiers while retaining monospace.
        if re.search(r"[/\\_.]", value):
            return protect(r"\path{" + value + "}")
        return protect(r"\texttt{\detokenize{" + value + "}}")

    text = re.sub(r"`([^`]+)`", render_code, text)
    text = re.sub(r"\*\*([^*]+)\*\*", lambda match: protect(r"\textbf{" + escape_plain(match.group(1)) + "}"), text)
    text = re.sub(r"\*([^*]+)\*", lambda match: protect(r"\emph{" + escape_plain(match.group(1)) + "}"), text)
    text = escape_plain(text)
    for index, value in enumerate(tokens):
        text = text.replace(escape_plain(f"@@TOKEN{index:04d}@@"), value)
    return text


def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def table_to_latex(lines: list[str], caption: str) -> list[str]:
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in lines]
    rows = [row for row in rows if not all(re.fullmatch(r":?-{3,}:?", cell) for cell in row)]
    column_count = len(rows[0])
    # Leave room for inter-column padding and use a narrow identifier column
    # where appropriate. Ragged-right cells avoid the extreme word spacing that
    # fully justified prose creates in narrow evidence tables.
    first_header = rows[0][0].lower()
    if column_count > 2 and first_header == "id":
        widths = [0.05] + [0.80 / (column_count - 1)] * (column_count - 1)
    elif column_count > 2 and first_header == "phase":
        widths = [0.09] + [0.76 / (column_count - 1)] * (column_count - 1)
    elif column_count == 2 and first_header == "artifact":
        widths = [0.52, 0.36]
    elif column_count == 2 and first_header == "phase":
        widths = [0.12, 0.76]
    else:
        widths = [0.88 / column_count] * column_count
    spec = "@{}" + " ".join(
        f">{{\\raggedright\\arraybackslash}}p{{{width:.3f}\\textwidth}}"
        for width in widths
    ) + "@{}"
    output = [r"\begin{longtable}{" + spec + "}", r"\caption{" + inline(caption) + r"}\\", r"\toprule"]
    output.append(" & ".join(r"\textbf{" + inline(cell) + "}" for cell in rows[0]) + r" \\")
    output.extend([r"\midrule", r"\endfirsthead", r"\toprule"])
    output.append(" & ".join(r"\textbf{" + inline(cell) + "}" for cell in rows[0]) + r" \\")
    output.extend([r"\midrule", r"\endhead"])
    for row in rows[1:]:
        rendered = []
        for cell in row:
            value = inline(cell)
            # File-system paths need discretionary breaks at slashes and dots;
            # \path is provided by hyperref/url and preserves monospace styling.
            value = re.sub(
                r"\\texttt\{\\detokenize\{([^{}]*[/\\\\][^{}]*)\}\}",
                lambda match: r"\\path{" + match.group(1) + "}",
                value,
            )
            rendered.append(value)
        output.append(" & ".join(rendered) + r" \\")
    output.extend([r"\bottomrule", r"\end{longtable}"])
    return output


def convert(lines: list[str], context: str) -> str:
    output: list[str] = []
    paragraph: list[str] = []
    table_number = 0
    current_heading = "Evidence summary"

    def flush() -> None:
        if paragraph:
            output.append(inline(" ".join(item.strip() for item in paragraph)))
            output.append("")
            paragraph.clear()

    index = 0
    in_math = False
    while index < len(lines):
        line = lines[index].rstrip()
        stripped = line.strip()
        if stripped == r"\[":
            flush(); in_math = True; output.append(r"\["); index += 1; continue
        if in_math:
            output.append(line)
            if stripped == r"\]":
                in_math = False; output.append("")
            index += 1; continue
        if not stripped:
            flush(); index += 1; continue
        heading = re.match(r"^(#{1,4})\s+(.+)$", stripped)
        if heading:
            flush()
            level, title = len(heading.group(1)), re.sub(r"^\d+(?:\.\d+)*\.?\s+", "", heading.group(2))
            current_heading = title
            if context == "appendices" and level == 2:
                title = re.sub(r"^Appendix\s+[A-Z]\.\s*", "", title)
                output.extend([r"\chapter{" + inline(title) + "}", r"\label{app:" + slug(title) + "}", ""])
            elif level == 1:
                output.extend([r"\chapter{" + inline(title) + "}", r"\label{chap:" + slug(title) + "}", ""])
            elif level == 2:
                output.extend([r"\section{" + inline(title) + "}", r"\label{sec:" + slug(title) + "}", ""])
            elif level == 3:
                output.extend([r"\subsection{" + inline(title) + "}", ""])
            else:
                output.extend([r"\subsubsection{" + inline(title) + "}", ""])
            index += 1; continue
        image_match = re.match(r"!\[([^\]]*)\]\(([^)]+)\)", stripped)
        if image_match:
            flush()
            caption = image_match.group(1)
            look = index + 1
            while look < len(lines) and not lines[look].strip():
                look += 1
            if look < len(lines) and lines[look].strip().startswith("**Figure"):
                caption = re.sub(r"^\*\*Figure\s+\d+\.\*\*\s*", "", lines[look].strip())
                index = look
            filename = Path(image_match.group(2)).name
            label = "fig:" + slug(Path(filename).stem)
            output.extend([
                r"\begin{figure}[p]", r"\centering",
                r"\includegraphics[width=0.96\textwidth,height=0.78\textheight,keepaspectratio]{figures/" + filename + "}",
                r"\caption{" + inline(caption) + "}", r"\label{" + label + "}", r"\end{figure}", "",
            ])
            index += 1; continue
        if stripped.startswith("|") and "|" in stripped[1:]:
            flush(); table_lines = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                table_lines.append(lines[index]); index += 1
            table_number += 1
            output.extend(table_to_latex(table_lines, f"{current_heading}: structured summary")); output.append("")
            continue
        if re.match(r"^[-*]\s+", stripped):
            flush(); items = []
            while index < len(lines) and re.match(r"^[-*]\s+", lines[index].strip()):
                items.append(re.sub(r"^[-*]\s+", "", lines[index].strip())); index += 1
            output.append(r"\begin{itemize}")
            output.extend(r"\item " + inline(item) for item in items)
            output.extend([r"\end{itemize}", ""]); continue
        if re.match(r"^\d+\.\s+", stripped):
            flush(); items = []
            while index < len(lines) and re.match(r"^\d+\.\s+", lines[index].strip()):
                items.append(re.sub(r"^\d+\.\s+", "", lines[index].strip())); index += 1
            output.append(r"\begin{enumerate}")
            output.extend(r"\item " + inline(item) for item in items)
            output.extend([r"\end{enumerate}", ""]); continue
        if stripped.startswith("```"):
            flush(); language = stripped[3:]; code = []; index += 1
            while index < len(lines) and not lines[index].strip().startswith("```"):
                code.append(lines[index]); index += 1
            output.extend([r"\begin{verbatim}", *code, r"\end{verbatim}", ""])
            index += 1; continue
        paragraph.append(stripped)
        index += 1
    flush()
    return "\n".join(output).rstrip() + "\n"


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    text = re.sub(r"\A---\n.*?\n---\n", "", text, flags=re.S)
    blocks = re.split(r"(?m)^#\s+", text)
    sections: dict[str, list[str]] = {}
    for block in blocks:
        if not block.strip():
            continue
        first, *rest = block.splitlines()
        sections[first.strip()] = ["# " + first.strip(), *rest]

    front_map = {
        "Declaration": "declaration.tex",
        "Acknowledgements": "acknowledgements.tex",
        "Abstract": "abstract.tex",
        "Abbreviations": "abbreviations.tex",
    }
    (TARGET / "frontmatter").mkdir(parents=True, exist_ok=True)
    for title, filename in front_map.items():
        content = convert(sections[title], "frontmatter")
        content = content.replace(
            r"\chapter{" + title + "}",
            r"\chapter*{" + title + r"}\addcontentsline{toc}{chapter}{" + title + "}",
            1,
        )
        (TARGET / "frontmatter" / filename).write_text(content, encoding="utf-8")

    (TARGET / "chapters").mkdir(parents=True, exist_ok=True)
    chapter_titles = [title for title in sections if re.match(r"^(?:[1-9]|10)\.\s", title)]
    for title in chapter_titles:
        number = int(title.split(".", 1)[0])
        clean = title.split(".", 1)[1].strip()
        filename = f"{number:02d}_{slug(clean)}.tex"
        (TARGET / "chapters" / filename).write_text(convert(sections[title], "chapter"), encoding="utf-8")

    (TARGET / "appendices").mkdir(parents=True, exist_ok=True)
    (TARGET / "appendices" / "appendices.tex").write_text(
        convert(sections["12. Appendices"][1:], "appendices"), encoding="utf-8"
    )
    print(f"Generated {len(chapter_titles)} chapters, {len(front_map)} front-matter files and appendices")


if __name__ == "__main__":
    main()
