-- PDF-only presentation: vector figures, portable links and joined captions.
local paper_dir = pandoc.path.directory(PANDOC_STATE.input_files[1])
-- The Markdown carries a readable H1; the PDF template already prints metadata title.
function Pandoc(doc)
  if doc.meta.title and #doc.blocks > 0 and doc.blocks[1].t == "Header"
     and doc.blocks[1].level == 1
     and pandoc.utils.stringify(doc.blocks[1].content) == pandoc.utils.stringify(doc.meta.title) then
    doc.blocks:remove(1)
  end
  return doc
end
local function project_path(path)
  local parts = {"docs", "paper"}
  for part in path:gmatch("[^/]+") do
    if part == ".." then
      assert(#parts > 0, "Reference escapes repository: " .. path)
      table.remove(parts)
    elseif part ~= "." then
      table.insert(parts, part)
    end
  end
  return table.concat(parts, "/")
end

function Image(el)
  local vector = el.src:gsub("%.png$", ".svg")
  local file = io.open(pandoc.path.join({paper_dir, vector}), "rb")
  if file then file:close(); el.src = vector end
  el.src = "/" .. project_path(el.src)
  el.attributes["width"] = "100%"
  return el
end

function Link(el)
  if not el.target:match("^%a[%w+.-]*:") and not el.target:match("^#") then
    local kind = el.target:match("/$") and "tree" or "blob"
    el.target = "https://github.com/peijin94/vamp/" .. kind .. "/main/" .. project_path(el.target)
  end
  return el
end

function Blocks(blocks)
  local result = pandoc.List()
  local i = 1
  while i <= #blocks do
    local item, caption = blocks[i], blocks[i+1]
    if item.t == "Para" and #item.content == 1 and item.content[1].t == "Image"
        and caption and caption.t == "Para" and pandoc.utils.stringify(caption):match("^Figure %d") then
      result:insert(pandoc.Div({item, caption}, pandoc.Attr("", {}, {["typst:breakable"]="false"})))
      i = i + 2
    else
      result:insert(item)
      i = i + 1
    end
  end
  return result
end

function Pandoc(doc)
  if doc.blocks[1].t == "Header" and doc.blocks[1].level == 1 then
    doc.meta.title = pandoc.MetaInlines(doc.blocks[1].content)
    doc.blocks:remove(1)
  end
  return doc:walk({Header = function(el) el.level = math.max(1, el.level-1); return el end})
end
