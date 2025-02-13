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
  -- relative fraction sizes for columns
  local fractions = {0.1429, 0.1905, 0.3333, 0.1429, 0.1905}

  for i, colspec in ipairs(tbl.colspecs) do
    -- Force alignment to left, set fraction
    colspec[1] = "AlignLeft"  -- Could be "AlignDefault", "AlignCenter", etc.
    colspec[2] = fractions[i]
  end

  return {
    pandoc.RawBlock("latex", "\\begin{landscape}"),
    tbl,
    pandoc.RawBlock("latex", "\\end{landscape}"),
  }
end
