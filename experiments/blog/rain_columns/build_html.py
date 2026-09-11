"""Insert the code excerpt from snippet.jl into blog.md and render blog.html from blog.md.

Usage: python3 build_html.py   (run from this directory)
The excerpt is the model definitions and the run calls; the full script is linked from the post.
"""
import html, re, pathlib

here = pathlib.Path(__file__).parent
snippet = (here / "snippet.jl").read_text().splitlines()

def block(start, stop):
    i = next(k for k, l in enumerate(snippet) if l.startswith(start))
    j = next(k for k in range(i + 1, len(snippet)) if snippet[k].startswith(stop))
    return snippet[i:j]

def strip_blank(b):
    while b and not b[-1].strip():
        b.pop()
    return b
excerpt = strip_blank(block("# --- the soil", "# --- initial state")) + [""] + strip_blank(block("# --- run both", "# --- plot"))
code = "\n".join(excerpt)

md = (here / "blog.md").read_text()
md = re.sub(r"<!-- SNIPPET -->\n```julia\n.*?```\n<!-- /SNIPPET -->",
            "<!-- SNIPPET -->\n```julia\n" + code + "\n```\n<!-- /SNIPPET -->", md, flags=re.S)
(here / "blog.md").write_text(md)

# ---------------------------------------------------------------- markdown → html (the subset used in blog.md)
def inline(t):
    t = html.escape(t, quote=False)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<em>\1</em>", t)
    t = re.sub(r"`(.+?)`", r"<code>\1</code>", t)
    t = re.sub(r"\[(.+?)\]\((.+?)\)", r'<a href="\2">\1</a>', t)
    return t

body = re.sub(r"<!--.*?-->\n?", "", md, flags=re.S)
lines = body.split("\n")
out, i, header = [], 0, {}
def wrap_open():
    out.append('<div class="wrap">')
def wrap_close():
    if out and out[-1] != "</div>": out.append("</div>")
wrap_open()
while i < len(lines):
    l = lines[i]
    if not l.strip():
        i += 1; continue
    if l.startswith("*CliMA software stack"):
        header["eyebrow"] = l.strip("*"); i += 1; continue
    if l.startswith("# "):
        header["h1"] = l[2:]; i += 1
        while not lines[i].strip(): i += 1
        header["dek"] = lines[i]; i += 1
        while not lines[i].strip(): i += 1
        header["byline"] = lines[i].strip("*"); i += 1
        out.append("<header>\n  <div class=\"eyebrow\">%s</div>\n  <h1>%s</h1>\n  <p class=\"dek\">%s</p>\n  <div class=\"byline\">%s</div>\n</header>"
                   % (inline(header["eyebrow"]), inline(header["h1"]), inline(header["dek"]), inline(header["byline"])))
        continue
    if l.startswith("## "):
        out.append("<h2>%s</h2>" % inline(l[3:])); i += 1; continue
    if l.startswith("!["):
        m = re.match(r"!\[(.*?)\]\((.*?)\)", l); i += 1
        while not lines[i].strip(): i += 1
        cap = lines[i].strip(); i += 1
        cap = re.sub(r"^\*\*\*(.+?)\*\*(.*)\*$", r"**\1**\2", cap)   # bold lead-in, then plain caption text
        cap = inline(cap).replace("<strong>", "<b>", 1).replace("</strong>", "</b>", 1)
        wrap_close()
        out.append('<div class="wide">\n<figure>\n  <img src="%s" alt="%s">\n  <figcaption>%s</figcaption>\n</figure>\n</div>'
                   % (m.group(2), html.escape(m.group(1), quote=True), cap))
        wrap_open(); continue
    if l.startswith("```"):
        j = i + 1
        while not lines[j].startswith("```"): j += 1
        out.append("<pre><code>%s</code></pre>" % html.escape("\n".join(lines[i + 1:j]), quote=False)); i = j + 1; continue
    if l.startswith("> "):
        paras = []
        while i < len(lines) and lines[i].startswith(">"):
            t = lines[i][1:].strip()
            if t: paras.append(t)
            i += 1
        title = paras.pop(0).strip("*")
        out.append('<div class="note">\n  <span class="eyebrow">%s</span>\n%s\n</div>' % (inline(title), "\n".join("  <p>%s</p>" % inline(p) for p in paras)))
        continue
    if l.strip() == "---":
        out.append("<hr>"); i += 1
        while not lines[i].strip(): i += 1
        out.append('<p class="foot">%s</p>' % inline(lines[i])); i += 1; continue
    out.append("<p>%s</p>" % inline(l)); i += 1
wrap_close()

css = (here / "style.css").read_text()
page = "<title>%s</title>\n%s\n\n%s\n" % (html.escape(header["h1"]), css, "\n\n".join(out))
(here / "blog.html").write_text(page)
print("wrote blog.md (snippet inserted) and blog.html")
