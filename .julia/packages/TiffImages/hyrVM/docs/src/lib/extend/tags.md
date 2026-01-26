# Built-in Tags

Tags are stored as an enum called `TiffTag`

```@docs
TiffImages.TiffTag
```

### Full list of built-in tags

```@eval
using TiffImages, Markdown
tags = instances(TiffImages.TiffTag)
mapping = collect.(zip(Int.(tags), string.(tags)))
insert!(mapping, 1, ["Tag ID", "Tag Description"])
Markdown.MD(Markdown.Table(mapping, fill(:l, length(mapping))))
```