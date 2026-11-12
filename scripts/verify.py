from pathlib import Path
from collections import Counter
import re,json,hashlib
R=Path(__file__).resolve().parents[1];seen={}
def visit(n):
 p=R/n
 if n in seen:return
 t=re.sub(r'(?<!\\)%[^\n]*','',p.read_text());seen[n]=t
 for q in re.findall(r'\\input\{([^}]+)\}',t):visit(q if q.endswith('.tex') else q+'.tex')
visit('papers/pfaffian_brill_noether_technical.tex')
t='\n'.join(seen.values());labels=re.findall(r'\\label\{([^}]+)\}',t);refs=re.findall(r'\\(?:ref|eqref|autoref|pageref)\{([^}]+)\}',t);bib=re.findall(r'\\bibitem(?:\[[^]]*\])?\{([^}]+)\}',t);cites=[x.strip() for g in re.findall(r'\\cite(?:\[[^]]*\])*\{([^}]+)\}',t) for x in g.split(',')]
assert len(labels)==len(set(labels));assert len(bib)==len(set(bib));assert not set(refs)-set(labels),set(refs)-set(labels);assert not set(cites)-set(bib),set(cites)-set(bib)
a={'tex_inputs':len(seen),'labels':len(labels),'references':len(refs),'citations':len(cites),'bibliography':len(bib),'sha256':{n:hashlib.sha256((R/n).read_bytes()).hexdigest() for n in seen},'limit':'Static closure only, not mathematical proof audit.'};(R/'build').mkdir(exist_ok=True);(R/'build/source_check.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2))
