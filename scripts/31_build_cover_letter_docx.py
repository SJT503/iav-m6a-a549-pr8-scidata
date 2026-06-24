# -*- coding: utf-8 -*-
"""31_build_cover_letter_docx.py — 构建 cover_letter.docx"""
import re, os
from pathlib import Path
import pypandoc

SRC = Path(r"D:\IAV_m6A\SciData_descriptor\submission\01_manuscript\cover_letter_source.md")
CLEAN = SRC.parent / "_cl_clean.md"
OUT = SRC.parent / "cover_letter.docx"

text = SRC.read_text(encoding="utf-8")
n = text.count("<!--")
text = re.sub(r"<!--.*?-->", "", text, flags=re.DOTALL)
text = re.sub(r"(?m)^\s*---\s*$", "", text)
text = re.sub(r"\n{3,}", "\n\n", text)
CLEAN.write_text(text, encoding="utf-8")
pypandoc.convert_file(str(CLEAN), "docx", outputfile=str(OUT),
    extra_args=["--from=markdown+pipe_tables"])
os.remove(str(CLEAN))
print(f"stripped {n} comments; saved {OUT.name} ({OUT.stat().st_size/1024:.1f} KB)")
