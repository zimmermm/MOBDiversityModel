@echo off
start /WAIT simstrat.exe json_rotsee_simstrat.par
echo run>> E:/polybox/SURF/MobTrait/results/diversity/instruction_server.txt
:check_instructions
sleep 1
(for /f "usebackq eol= " %%a in ("E:/polybox/SURF/MobTrait/results/diversity/instruction_server.txt") do break) && goto :check_instructions || echo done