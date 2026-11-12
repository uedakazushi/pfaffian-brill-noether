# Pfaffian–Brill–Noether境界・接触分類

2026-11-13 / v1。正本は本リポジトリのmain。

現行入口は `papers/pfaffian_brill_noether_technical.tex`。三次三様体のBrill–Noether境界、六次Pfaffian接錐・解消・Torelli、全次数接触商・頂点核、Tate/Hessian補正、同変代数化、五変数Hessian特異点を扱う。まず完全な技術稿を保持し、分割は依存監査後に行う。

```sh
latexmk papers/pfaffian_brill_noether_technical.tex
python3 scripts/verify.py
python3 computations/verify_pfaffian_vertex_and_tate.py
```

日本語LuaLaTeX・LuaTeX-ja・HaranoAjiフォント・latexmk、Python3とSymPyが必要。PDFと中間物はbuild/で非追跡。VS Codeはリポジトリをフォルダとして開き、既定のLaTeX Workshop recipeを使う。Actions・force pushなし。

一般Poisson分類・量子化はnc-cubic論文Aの既知入力であり、本稿の独立新成果と数えない。原稿内に証明があること、有限検算、独立監査、新規性を区別する。欠けた旧大型行列と未監査範囲は `notes/research_notes.md`、実行範囲は `validation.json`。
