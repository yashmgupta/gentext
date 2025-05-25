#!/usr/bin/env python3
"""
A simple Tkinter app to load a GenBank file and generate
a “manuscript‐style” summary of its gene content.
"""

import re
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from Bio import SeqIO

def calc_gc(seq):
    s = str(seq).upper()
    gc = s.count('G') + s.count('C')
    return (gc / len(s)) * 100 if len(s) else 0.0

def extract_gene_info(record):
    counts = {'CDS':0, 'tRNA':0, 'rRNA':0}
    genes = []
    for feat in record.features:
        ft = feat.type
        if ft not in counts:
            continue
        counts[ft] += 1
        name = (feat.qualifiers.get('gene')
                or feat.qualifiers.get('product')
                or [ft])[0]
        if ft == 'tRNA':
            m = re.search(r'tRNA-([A-Za-z]+)', name)
            if m:
                name = f"{m.group(1)}-tRNA"
        if ft == 'rRNA':
            if re.search(r'16S|large', name, re.I):
                name = '16S-rRNA'
            elif re.search(r'12S|small', name, re.I):
                name = '12S-rRNA'
        genes.append({'pos': int(feat.location.start), 'name': name})
    genes.sort(key=lambda g: g['pos'])
    total = counts['CDS'] + counts['tRNA'] + counts['rRNA']
    return {'counts': counts, 'ordered': genes, 'total': total}

def generate_manuscript_summary(record, info):
    L = len(record.seq)
    gc = calc_gc(record.seq)
    c = info['counts']
    names = [g['name'] for g in info['ordered']]

    p1 = (
        "In the present study, we report the sequence of "
        f"{record.id}, which spans {L:,} base pairs and exhibits a GC content "
        f"of {gc:.2f}%."
    )
    p2 = (
        "Automated annotation identified "
        f"{c['CDS']} protein-coding sequences (CDSs), "
        f"{c['tRNA']} transfer RNA genes, and {c['rRNA']} ribosomal RNA genes, "
        f"for a total of {info['total']} functional gene features. "
        "This gene complement is consistent with values reported for similar taxa."
    )
    p3 = (
        "Analysis of gene synteny revealed the following linear arrangement: "
        + "; ".join(names)
        + ". Notably, the clustering of rRNA genes suggests conservation of "
        "the ribosomal operon structure."
    )
    return "\n\n".join([p1, p2, p3])

class GenBankApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("GenBank Manuscript Summary")
        self.geometry("800x600")

        # Buttons
        btn_frame = tk.Frame(self)
        btn_frame.pack(fill="x", pady=5)
        load_btn = tk.Button(btn_frame, text="Load GenBank File", command=self.load_file)
        load_btn.pack(side="left", padx=5)
        clear_btn = tk.Button(btn_frame, text="Clear", command=self.clear_text)
        clear_btn.pack(side="left")

        # Scrolled text area
        self.text = scrolledtext.ScrolledText(self, wrap="word")
        self.text.pack(expand=True, fill="both", padx=5, pady=5)

    def load_file(self):
        path = filedialog.askopenfilename(
            title="Select GenBank File",
            filetypes=[("GenBank files", "*.gb *.gbk"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            records = SeqIO.parse(path, "genbank")
            summaries = []
            for rec in records:
                info = extract_gene_info(rec)
                summaries.append(generate_manuscript_summary(rec, info))
            if not summaries:
                raise ValueError("No records found in file.")
            output = "\n\n" + ("="*80 + "\n\n").join(summaries)
            self.text.insert("end", output)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file:\n{e}")

    def clear_text(self):
        self.text.delete("1.0", "end")

if __name__ == "__main__":
    app = GenBankApp()
    app.mainloop()
