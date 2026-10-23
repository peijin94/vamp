#set document(title: "$title$", author: "", keywords: ("Measurement Set", "visibility compression", "quantization"))
#set page(paper: "a4", margin: (x: 22mm, top: 20mm, bottom: 20mm), numbering: "1")
#set text(font: "DejaVu Serif", size: 10pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: none)
#show heading: set par(justify: false)
#show heading: set text(hyphenate: false)
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 11.5pt, weight: "bold")
#show link: set text(fill: rgb("#23618d"))
#show raw: set text(font: "DejaVu Sans Mono", size: 8pt)
#show table: set text(size: 8.5pt)
#show table: set par(justify: false)
#set table(inset: 4pt, stroke: 0.3pt + rgb("#d8dee4"))
#show table.cell.where(y: 0): set text(weight: "bold")
#let horizontalrule = line(length: 100%, stroke: 0.5pt + rgb("#d8dee4"))
#align(center)[
  #set par(justify: false)
  #text(size: 18pt, weight: "bold", hyphenate: false)[$title$]
]
#v(0.6em)
$body$
