
-- string functions not yet availible in Foldit
upper = {a="A", c="C", d="D", e="E", f="F", g="G", h="H", i="I", k="K", l="L", m="M", n="N", p="P", q="Q", r="R", s="S", t="T", v="V", w="W", x="X", y="Y"}

function wiggle_side_all()
	do_global_wiggle_sidechains(5)
	do_global_wiggle_all(1)
end

function wiggle_sphere(radius, resi)
	while 1 do
		wiggle_side_all()
        	if confirm("Keep wiggling?") then i=1 else break end
	end
end

function mutation_trial (original, present_aa)
	reset_recent_best()
        before=get_score()
        deselect_all()					-- Redo deselect/select here
        select_segment_pdb_index("A", original[2])	-- As user may have changed things in examination
        replace_aa(original[1])
        wiggle_sphere(10, original[2])
        after=get_score()
        if confirm("Keep the change? ("..present_aa..original[2]..original[1]..")\nScore delta is: \n"..after-before.."\nNegative score means native is favored. Click yes to keep the change.") then keep="yes" else keep="no" end
	if keep=="no" then restore_recent_best() end
end

function ask_mutation(original, present_aa)
	if confirm("Try this mutation at "..original[2].."? \nCurrent residue is "..present_aa.."\nOriginal residue is "..original[1]) then mutation_trial(original, present_aa) end
end

for key, original in pairs(parent_table) do 
        deselect_all()
        select_segment_pdb_index("A", original[2])
    	present_aa = upper[get_aa(original[2])]
	if present_aa ~= original[1] then ask_mutation(original, present_aa) end
end

deselect_all() -- Remove selection from last amino acid
