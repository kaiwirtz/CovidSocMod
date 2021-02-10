# Parameter file written by XSiSi 2.0
# Last modified at Tue Feb 6 11:11:11 CET 2021
comment filenames for output
string	MovieResFile	"res/ResRef_0.outb"
	\d file name for binary simulation results
comment Control parameter of variation
string	VarResFile	"res/No_00.var"
	\d name of ASCII file for results of parameter variation ForceLockTime
string	VarIndexFile	"varindex1a.bin"
	\d name of binary file for reading best indices of Parameter Variation
comment Control parameters for experiments (II)
int	VarActive	0
	\d 1: parameter variation on/off
int	NumDice		5920
	\d 10: output during variation <0: emulates random parameter set for check of variations
int	MonteCarlo	0
	\d 0: simultaneous var 1: Random parameter variation 2: OAT (single) variation
int	VarOutputStep	0
	\d PDF every 2^n output, <0 no PD
int	AnalysisType	11
	\d 0:
array parindex
typeOfArray int
dimension 2 2
 data
        4 0
        10 9365
end
comment Data integration and comparison
int     DataActive      1
        \d 1: comparison with data
int	DataCountryIndex	3
				\d 0: all together; <0 loop over countries; >0 single country
comment 1:Belgium 2:France 3:Germany 4:Iran 5:Ireland 6:Italy 7:Netherlands 8:Portugal 9:Spain 10:Sweden 11:Switzerland 12:United Kingdom 13:Georgia 14:Indiana 15:Louisiana 16:Massachusetts 17:Michigan 18:New Jersey 19:New York 20:Washington
float	critfatal	6E-7
   \d critical fatality when infection is restarted in simulation and data compared
array   err_data_weights
        \d 0:fatality
        typeOfArray     float
        dimension       1
data
        1
end
float   data_off_tim    0
        \u days
        \d SimTime without data comparison  -34.
comment Control parameters for experiments (II)
int     Experiment   23
        \d <-1: loop over IndExp;
comment Scenarios
array   IndExp
        \d Index of experiment starting from 0
        typeOfArray     int
        dimension       9
data
	-1 10 60 155 250 270 366 -1 -1
end
comment output parameters
float	storetim	0.0
	\u days
	\d Sim time without output from TimeStart
comment pseudo zero for floating point arithmetic
float	Zero	1.0E-6
comment numerical settings
float	mindelt	1.0E-9
float	maxdelt	1.0
float	relrate	1.0
float	relchange	0.2
float	accuracy	0.1
int	method	1
	\d Integration type 0:Euler 1:RungeKutta2-3
