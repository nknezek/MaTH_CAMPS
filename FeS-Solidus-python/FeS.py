mol_Fe = 55.845;
mol_FeS = 87.91;

'''
FeS.py
Methods  to convert between sulfur weight percent and mol fraction in an Fe-FeS liquid
system.
By Nicholas Knezek supported by CIDER
July 2014
'''

def mf_to_wpS(FeSfrac):
	FeSwp = FeSfrac*mol_FeS/(FeSfrac*mol_FeS + (1.-FeSfrac)*mol_Fe);
	Swp = FeSwp*(mol_FeS-mol_Fe)/mol_FeS;
	return Swp

def wpS_to_mf(Swp):
	FeSwp = Swp*mol_FeS/(mol_FeS-mol_Fe)
	FeSfrac = FeSwp*mol_Fe/(mol_FeS + FeSwp*(mol_Fe-mol_FeS))
	return FeSfrac

###### Test Conversions #######
# S = 0.364748
# mfS = 1.0
# print wpS_to_mf(S)
# print mf_to_wpS(mfS)