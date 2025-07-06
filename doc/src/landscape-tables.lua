-- helper: wrap any block in \begin{landscape}..\end{landscape}
local function in_landscape(block)
  return {
    pandoc.RawBlock("latex", "\\begin{landscape}\\begingroup\\footnotesize"),
    block,
    pandoc.RawBlock("latex", "\\endgroup\\end{landscape}")
  }
end

----------------
--  NEW CODE  --
----------------
-- Catch Pandoc 3.x SimpleTable
function SimpleTable(tbl)
  if not FORMAT:match("latex")   then return tbl end
  if #tbl.align ~= 5             then return tbl end
  return in_landscape(tbl)
end

-- old pandoc table type
function Table(tbl)
  if not FORMAT:match("latex") then return tbl end

  return in_landscape(tbl)
end
