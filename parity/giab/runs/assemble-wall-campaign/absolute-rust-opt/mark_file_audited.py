#!/usr/bin/env python3
"""Usage: mark_file_audited.py <path> [verdict=tight] [note=...]
Sets all pending fns in path to verdict; marks AUDIT_BOARD audited when none pending.
"""
import sys
from pathlib import Path
root = Path(__file__).resolve().parent
path = sys.argv[1]
verdict = "tight"
note = "reviewed isolation"
for a in sys.argv[2:]:
    if a.startswith("verdict="): verdict = a.split("=",1)[1]
    if a.startswith("note="): note = a.split("=",1)[1]
vpath = root/"FN_VERDICTS.tsv"
lines = vpath.read_text().splitlines()
hdr, body = lines[0], lines[1:]
out=[]
for ln in body:
    parts=ln.split("\t")
    if parts[1]==path and parts[4]=="pending":
        parts[4]=verdict
        parts[5]=note
        out.append("\t".join(parts))
    else:
        out.append(ln)
vpath.write_text(hdr+"\n"+"\n".join(out)+"\n")
# board
board=root/"AUDIT_BOARD.tsv"
bl=board.read_text().splitlines()
bh, bb = bl[0], bl[1:]
pending=any(ln.split("\t")[1]==path and ln.split("\t")[4]=="pending" for ln in out)
nb=[]
for ln in bb:
    p=ln.split("\t")
    if p[1]==path:
        p[3]="audited" if not pending else p[3]
        nb.append("\t".join(p))
    else:
        nb.append(ln)
board.write_text(bh+"\n"+"\n".join(nb)+"\n")
print(path, "done pending=", pending)
