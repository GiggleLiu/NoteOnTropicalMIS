using TropicalNumbers, OMEinsum, NoteOnTropicalMIS

code = ein"u,v,w,j,k,vw,vj,vk,ku,wu,wk->juk"
mis_contract(Tropical(1.0), code)