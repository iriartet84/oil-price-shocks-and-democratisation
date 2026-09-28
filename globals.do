/*------------------------------------------------------------------
M16 pipeline — shared paths. `do`-ed at the top of every script.
------------------------------------------------------------------*/
* Every path is relative to the project folder (the one holding Data/ and
* Code/). Each script cd's there before calling this file, so the package
* runs unchanged on any machine, whether Stata starts in the project folder
* or in Code/ (e.g. after double-clicking a .do file).
global root     "`c(pwd)'"
global data     "$root/Data/data.dta"
global balanced "$root/Data/balanced_conditional_curse.dta"
global out_orig "$root/Output/original"
global out_ext  "$root/Output/extended"
global logs     "$root/Output/logs"

capture mkdir "$root/Output"
capture mkdir "$out_orig"
capture mkdir "$out_ext"
capture mkdir "$logs"

* xtabond2's cluster() option requires Mata's speed-favoring mode.
* Session-only (not `perm`) so we don't alter the shared network install's defaults.
mata: mata set matafavor speed
