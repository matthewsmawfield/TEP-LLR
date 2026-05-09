C  met value storage
	COMMON /MET/ met_sod, met_p, met_t, met_h, met_o, n_met, 
     &               n_met_next
	integer n_met, n_met_next, met_o(500)
        double precision met_sod(500), met_p(500), met_t(500),
     +          met_h(500)
