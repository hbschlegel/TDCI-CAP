-- colwidths-landscape.lua
-- Lua filter for Pandoc that:
--   1. Forces 5-column Markdown tables to specific column widths
--   2. Wraps tables in a landscape environment (LaTeX/PDF only)
--
-- Usage example:
--   pandoc yourfile.md -o out.pdf --pdf-engine=xelatex \
--       -V tables=yes -V longtable \
--       --lua-filter=colwidths-landscape.lua

function Table(tbl)
  -- Only modify if output is LaTeX (PDF)
  if not FORMAT:match("latex") then
    return tbl
  end

  -- If this table doesn't have exactly 5 columns, do nothing
  local ncols = #tbl.colspecs
  if ncols ~= 5 then
    return {
      pandoc.RawBlock("latex", "\\begin{landscape}"),
      tbl,
      pandoc.RawBlock("latex", "\\end{landscape}"),
    }
  end

  ----------------------------------------------------------------------------
  -- 1) Force each of the 5 columns to a chosen fraction of the total width
  --    Pandoc's colspecs field is an array of (Alignment, Width),
  --    where Width is a fraction from 0.0 to 1.0 of total table width.
  ----------------------------------------------------------------------------

  -- Suppose you want these approximate relative widths:
  --   Column 1 = 3   cm
  --   Column 2 = 4   cm
  --   Column 3 = 7   cm
  --   Column 4 = 3   cm
  --   Column 5 = 4   cm
  -- Total = 3+4+7+3+4 = 21 cm
  -- So the fractions are 3/21, 3/21, 6/21, 3/21, 6/21

  local fractions = {0.1429, 0.1429, 0.2857, 0.1429, 0.2857}

  for i, colspec in ipairs(tbl.colspecs) do
    -- Force alignment to left, set fraction
    colspec[1] = "AlignLeft"  -- Could be "AlignDefault", "AlignCenter", etc.
    colspec[2] = fractions[i]
  end

  ----------------------------------------------------------------------------
  -- 2) Wrap the table in a landscape environment
  ----------------------------------------------------------------------------

  return {
    pandoc.RawBlock("latex", "\\begin{landscape}"),
    tbl,
    pandoc.RawBlock("latex", "\\end{landscape}"),
  }
end
