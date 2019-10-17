
def write_out_mse_clima(imax, jmax, zmax,mse2, mse2_adv, mse2_div, mdiv2, madv2, tadv2,  omse2, prefixout):
	mse2 = mse2.T
	
	mse2_adv = mse2_adv.T
	
	mse2_div = mse2_div.T
	
	mdiv2 = mdiv2.T
	
	madv2 = madv2.T
	
	tadv2 = tadv2.T
	
	omse2 = omse2.T
	
##  write out full MSE variables 
	nameout = prefixout + "MSE_mse_clim.out"
	fh = open(nameout, "wb")
	mse2.tofile(fh)
	fh.close()

	nameout = prefixout + "MSE_adv_clim.out"
	fh = open(nameout, "wb")
	mse2_adv.tofile(fh)
	fh.close()

	nameout = prefixout + "MSE_div_clim.out"	
	fh = open(nameout, "wb")
	mse2_div.tofile(fh)
	fh.close()

	nameout = prefixout + "MSE_mdiv_clim.out"
	fh = open(nameout, "wb")
	mdiv2.tofile(fh)
	fh.close()

	nameout = prefixout + "MSE_madv_clim.out"
	fh = open(nameout, "wb")
	madv2.tofile(fh)
	fh.close()

	nameout = prefixout + "MSE_tadv_clim.out"
	fh = open(nameout, "wb")
	tadv2.tofile(fh)
	fh.close()

	nameout = prefixout + "MSE_omse_clim.out"
	fh = open(nameout, "wb")
	omse2.tofile(fh)
	fh.close()	

