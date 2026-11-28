"""Build the editable journal manuscript and human-readable repository source.

Requires pandoc with citeproc and python-docx. Does not read experimental data.
The PDF is rendered separately from the final DOCX for visual verification.
"""
from pathlib import Path
import argparse
from io import BytesIO
import json
import re
import subprocess
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

HERE=Path(__file__).resolve().parent

def set_font(style,size=11):
    style.font.name='Times New Roman'
    style.font.size=Pt(size)
    style.font.color.rgb=RGBColor(0,0,0)
    rpr=style.element.get_or_add_rPr()
    fonts=rpr.find(qn('w:rFonts'))
    if fonts is None: fonts=OxmlElement('w:rFonts');rpr.insert(0,fonts)
    for k in ('ascii','hAnsi','eastAsia','cs'):fonts.set(qn('w:'+k),'Times New Roman')

def reference_doc(path):
    # Preserve Pandoc's native table/Compact styles and document settings.
    # A fresh python-docx template lacks these and some Word-compatible
    # renderers misinterpret tables produced against that incomplete reference.
    d=Document(BytesIO(subprocess.check_output(
        ['pandoc','--print-default-data-file','reference.docx'])))
    # Pandoc capitalizes built-in heading names; python-docx's name aliases
    # assume the lower-case spelling. Resolve by the actual style objects.
    styles={st.name:st for st in d.styles}
    for name in ('Header','Footer'):
        if name not in styles:
            styles[name]=d.styles.add_style(name,WD_STYLE_TYPE.PARAGRAPH)
    s=d.sections[0]
    s.page_width=Inches(8.5);s.page_height=Inches(11)
    s.top_margin=Inches(.72);s.bottom_margin=Inches(.72)
    s.left_margin=Inches(.85);s.right_margin=Inches(.85)
    s.header_distance=Inches(.28);s.footer_distance=Inches(.28)
    for name in ('Normal','Body Text','First Paragraph'):
        if name not in styles:continue
        st=styles[name];set_font(st,11)
        st.paragraph_format.space_after=Pt(6)
        st.paragraph_format.line_spacing=1.13
        st.paragraph_format.widow_control=True
    for name,size in [('Title',18),('Subtitle',11),('Heading 1',14),('Heading 2',12),('Heading 3',11)]:
        st=styles[name];set_font(st,size)
        st.font.bold=name!='Subtitle'
        st.paragraph_format.space_before=Pt(12 if name.startswith('Heading') else 0)
        st.paragraph_format.space_after=Pt(6)
        st.paragraph_format.keep_with_next=True
    set_font(styles['Caption'],9.5)
    styles['Caption'].paragraph_format.line_spacing=1.05
    styles['Caption'].paragraph_format.space_after=Pt(8)
    for nm in ('Header','Footer'):
        set_font(styles[nm],9)
    # Page numbers provide useful navigation for a manuscript.
    p=s.footer.paragraphs[0];p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');p._p.append(field)
    d.save(path)

def shade(cell,fill):
    pr=cell._tc.get_or_add_tcPr();e=pr.find(qn('w:shd'))
    if e is None:e=OxmlElement('w:shd');pr.append(e)
    e.set(qn('w:fill'),fill)

def format_document(path):
    d=Document(path)
    for style in d.styles:
        if style.type==1 and (style.name.startswith('Heading') or style.name in ('Title','Subtitle','Header')):
            style.font.color.rgb=RGBColor(0,0,0)
    for p in d.paragraphs:
        p.paragraph_format.widow_control=True
        if p._p.xpath('.//w:drawing'):
            p.alignment=WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.keep_with_next=True
            p.paragraph_format.space_before=Pt(5)
        if p.style.name in ('Caption','Image Caption') or p.text.startswith('Figure '):
            p.paragraph_format.keep_together=True
            p.paragraph_format.space_after=Pt(9)
            for r in p.runs:r.font.size=Pt(9.5)
        if p.text.startswith(('Table ','Supplementary Table ')):
            p.paragraph_format.keep_with_next=True
            p.paragraph_format.space_before=Pt(7)
            for r in p.runs:r.font.size=Pt(9.5)
        if p.style.name.startswith('Heading') or p.style.name=='Title':
            for r in p.runs:r.font.color.rgb=RGBColor(0,0,0)
            pr=p._p.get_or_add_pPr()
            for b in list(pr.findall(qn('w:pBdr'))):pr.remove(b)
        if p.style.name=='Bibliography':
            p.paragraph_format.space_after=Pt(6)
            p.paragraph_format.line_spacing=1.05
            for r in p.runs:r.font.size=Pt(9.5)
    for t in d.tables:
        t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
        n=len(t.columns)
        # Named features and protocol history need a larger prose column.
        if n==2: widths=[1.25,5.55]
        elif n==3: widths=[3.1,1.6,2.1]
        elif n==4: widths=[2.05,1.4,1.55,1.8]
        elif n==5: widths=[2.45,1.09,1.09,1.09,1.08]
        else: widths=[6.8/n]*n
        header=' '.join(c.text for c in t.rows[0].cells)
        if 'Primary AP estimand' in header:widths=[1.35,2.15,.75,1.65,.9]
        if 'Protocol record' in header:widths=[1.35,1.25,1.7,2.5]
        if 'TF-varying' in header:widths=[3.4,1.7,1.7]
        for i,col in enumerate(t.columns):col.width=Inches(widths[i])
        pr=t._tbl.tblPr
        borders=pr.find(qn('w:tblBorders'))
        if borders is None:borders=OxmlElement('w:tblBorders');pr.append(borders)
        for k in ('top','left','bottom','right','insideH','insideV'):
            e=OxmlElement('w:'+k);e.set(qn('w:val'),'single');e.set(qn('w:sz'),'4');e.set(qn('w:color'),'D9D9D9');borders.append(e)
        for ri,row in enumerate(t.rows):
            rpr=row._tr.get_or_add_trPr()
            e=OxmlElement('w:cantSplit');rpr.append(e)
            if ri==0:
                repeat=OxmlElement('w:tblHeader');rpr.append(repeat)
            for ci,cell in enumerate(row.cells):
                cell.width=Inches(widths[ci]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
                shade(cell,'263A49' if ri==0 else ('F2F5F7' if ri%2==0 else 'FFFFFF'))
                cp=cell._tc.get_or_add_tcPr();mar=OxmlElement('w:tcMar')
                for k,v in [('top','85'),('bottom','85'),('left','90'),('right','90')]:
                    el=OxmlElement('w:'+k);el.set(qn('w:w'),v);el.set(qn('w:type'),'dxa');mar.append(el)
                cp.append(mar)
                for p in cell.paragraphs:
                    p.paragraph_format.space_after=Pt(0);p.paragraph_format.space_before=Pt(0)
                    p.paragraph_format.line_spacing=1.05
                    p.paragraph_format.keep_with_next=False
                    p.alignment=WD_ALIGN_PARAGRAPH.LEFT if ci==0 or n==2 or 'Protocol record' in header or 'Primary AP estimand' in header else WD_ALIGN_PARAGRAPH.CENTER
                    for r in p.runs:
                        r.font.name='Times New Roman';r.font.size=Pt(9.2)
                        r.font.color.rgb=RGBColor(255,255,255) if ri==0 else RGBColor(0,0,0)
                        if ri==0:r.font.bold=True
    d.core_properties.title='A Controlled Evaluation of Graph Conditioned Evidence Weighting for Gene Regulatory Network Inference'
    d.core_properties.subject='Controlled GCEV evaluation and source-use diagnostics'
    d.core_properties.author=''
    d.core_properties.comments=''
    d.save(path)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    source=HERE/'manuscript_source.md'
    text=source.read_text()
    for marker in ('INSERTION_POINT','TODO','[M1','[S1]'):
        if marker in text:raise ValueError(f'Unresolved marker: {marker}')
    keys=set(re.findall(r'@([A-Za-z][A-Za-z0-9_]*)',text))
    refs=json.loads((HERE/'references.csl.json').read_text())
    missing=keys-{r['id'] for r in refs}
    if missing:raise ValueError(f'Missing references: {missing}')
    ref=out/'reference_style.docx';reference_doc(ref)
    common=['pandoc',str(source),'--from=markdown','--standalone','--citeproc',
            '--bibliography='+str(HERE/'references.csl.json'),'--csl='+str(HERE/'numeric.csl'),
            '--resource-path='+str(HERE)]
    target=out/'GCEV_Journal_Manuscript_Draft_v0_1_2026_11_28.docx'
    subprocess.run(common+['--reference-doc='+str(ref),'-o',str(target)],check=True)
    format_document(target)
    subprocess.run(common+['--to=gfm','--wrap=none','-o',str(HERE/'MANUSCRIPT.md')],check=True)
    # A self-contained LaTeX manuscript preserves editable scientific equations.
    subprocess.run(common+['--to=latex','-V','geometry:margin=0.85in','-V','fontsize=11pt',
                          '-V','colorlinks=true','-V','linkcolor=black','-V','urlcolor=blue',
                          '-o',str(HERE/'manuscript.tex')],check=True)
    print(json.dumps({'docx':str(target),'references_cited':len(keys),'tables':len(Document(target).tables),
                      'inline_figures':len(Document(target).inline_shapes),'source_words':len(text.split())}))

if __name__=='__main__':main()
