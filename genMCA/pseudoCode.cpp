/*procedure Inner_Differential_Evolution
begin
	intiate(M)
	repeat
		mutation operation
		crossover operation
		selesction operation,get newM
		M <- newM
	until stopcriterion
	find the best agent b in M
	r' <- b
	Return r'
end


procedure Outer_Differential_Evolution(degree m)
begin
	intiate(M')
	repeat
		mutation operation
		crossover operation
		selesction operation
			repeat
				m rows <- m groups <- each agent
				foreach i=1 to m
					find the best ri through Inner_DE
			until all agents have been dealt with
			get newM'
		M' <- newM'
	until stopcriterion
	find the best agent b in M'
	Return distribution b
end

procedure Nested_Differential_Evolution
begin
	repeat
		select one row r from MCA
		set initial solution s0 to r
		calculate ur, T(r)
		set degree d=1
		repeat
			find d rows to replace r
				use Outer_DE(d),get a distribution b
				m rows <- b
			iff T(s0)>T(m rows)
				s0 <- m rows
			d <- d+1
		until d>ur
		update MCA, replace r with s0
	until can not find one row that can be replaced
end


*/