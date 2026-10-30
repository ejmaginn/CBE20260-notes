# CBE 20260 Course Notes

Lecture notes and demo notebooks for CBE 20260, Introduction to Chemical
Engineering Thermodynamics. The rendered site is at
<https://ejmaginn.github.io/CBE20260-notes/>.

## Software

Everything except TeX comes from the pixi environment defined in
`pixi.toml` (Python, numpy, scipy, pandas, matplotlib, CoolProp,
ipywidgets, ipympl, JupyterLab, Quarto). PDF output uses the system TeX
installation (MacTeX).

First time on a new machine:

```bash
gh repo clone ejmaginn/CBE20260-notes Lecture-Notes-2026
cd Lecture-Notes-2026
pixi install
```

## Everyday commands

Run these from this folder.

| Command | What it does |
|---|---|
| `pixi run preview` | Live preview in the browser; updates as you save edits |
| `pixi run render-html` | Build the web pages only, into `_site/` (fast) |
| `pixi run render` | Build the web pages and the PDFs |
| `pixi run lab` | Open JupyterLab to work on the demo notebooks |
| `pixi run quarto render class09/lecture09.qmd --to pdf` | Build one lecture's PDF |
| `pixi run quarto publish gh-pages` | Render and publish the site to GitHub Pages |
| `pixi shell` | Open a shell with the environment active (`exit` to leave) |

To add a Python package, run `pixi add <package>` and commit the
changed `pixi.toml` and `pixi.lock`.

## Working from two laptops

GitHub is the hand-off point. Start a session with `git pull` (and
`pixi install`, which only does work if `pixi.lock` changed); end it with
`git add`, `git commit` and `git push`.

## Layout

- `classNN/`: one folder per class, with the lecture `.qmd`, its images,
  and any demo notebooks. Notebooks are rendered from their saved
  outputs, so re-run and save a notebook to update what the site shows.
- `index.qmd`: the site's home page.
- `_quarto.yml`: site settings, navbar and sidebar entries, PDF format.
- `latex/preamble.tex`: LaTeX macros for the PDF output.
- `references.bib`: bibliography.
