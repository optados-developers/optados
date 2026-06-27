#!/usr/bin/env python3

"""
J. A. J. Whaley-Baldwin, June 2026

Plotting routines for various OptaDOS vib-EELS outputs.
"""

import argparse
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm,LinearSegmentedColormap,PowerNorm
from matplotlib.ticker import AutoMinorLocator, NullLocator
from scipy.ndimage import gaussian_filter1d

##########################################################################################
# Global mpl plot settings.

plt.rcParams['axes.formatter.useoffset'] = False
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

##########################################################################################
# Utility functions.

def clean_line(line):
	line = line.replace("\n","")
	line = [ x for x in line.split(" ") if x ]
	return line

##########################################################################################
# Impact vib-EELS heatmap, optionally with vertical slices at certain q-points.

def plot_impact_heatmap(seedname, hsp_labels=None, qslice_fracs=None, elim=None):

	if (  (qslice_fracs) and (len(qslice_fracs) > 4) ):
		print("ERROR: Maximum number of q-slices for impact_heatmap is 4")
		exit()

	file = open("./"+seedname+"_impact-heatmap.dat","r")
	lines = file.readlines()

	energies = []
	read_data = False
	read_hsp_locs = False
	for line in lines:
		line = clean_line(line)
		# Get N_qpts.
		if ( ("Number" in line) and ("q-points" in line) ):
			N_qpts = int(line[-1])
			continue
		# Get N_energies.
		if ( ("Number" in line) and ("energies" in line) ):
			N_energies = int(line[-1])
			continue
		if ("High-symmetry" in line):
			read_hsp_locs = True
			continue
		if (read_hsp_locs):
			hsp_locs = [ int(x) - 1 for x in line ]
			read_hsp_locs = False
			continue
		if ("range" in line):
			continue
		# Now, allocate space for heatmap data, and parse.
		if ( ("Total" in line) and ("Loss" in line) ):
			heatmap_data = np.zeros( (N_energies,N_qpts) )
			qpt_idx_old = -1
			energy_count = 0
			read_data = True
			continue
		if (read_data):
			if ( len(line) == 0 ):
				continue
			else:
				qpt_idx = int(line[0])
				energy = float(line[1])
				loss = float(line[2])
				if (qpt_idx != qpt_idx_old):
					energy_count = 0
				heatmap_data[energy_count][qpt_idx-1] = loss
				if ( len(energies) < N_energies ):
					energies.append(energy)
				energy_count += 1
				qpt_idx_old = qpt_idx
				if ( (qpt_idx == N_qpts) and (energy_count == N_energies) ):
					break

	energies = np.array(energies)
	max_energy = np.max(energies)

	# Now plot the heatmap, and if requested, the q-slices.

	# If doing q-slices, create plot mosaic.
	if (qslice_fracs):
		qslice_colors = ["red","magenta","darkorange","forestgreen",]
		slices_todo = [ int(x * N_qpts) for x in qslice_fracs ]
		N_slices = len(slices_todo)
		mosaic = [ ["left", f"right_{i}"] for i in range(N_slices) ]
		# The total height of the figure changes to accommodate the number of requested q-slices.
		# So, set sensible min/max heights.
		fig_height = max(6.0, 2.2*N_slices)
		fig_height = min(7.0, fig_height)
		fig,ax = plt.subplot_mosaic( mosaic,figsize=(10,fig_height),layout="constrained",gridspec_kw={"width_ratios":[2,1],"wspace": 0.05} )
	# Otherwise, just create a single ax.
	else:
		fig,ax = plt.subplots(1,1)
		ax = {"left":ax}

	# Optionally blur in the x-direction, to simulate finite q-resolution.
	# Here, 'sigma_x' is in units of x-axis pixels, so it corresponds to your data step size.
	# Very roughly, for a path spacing of 0.01, a sigma_x of around 15 is sensible for most systems.
	#sigma_x = 15
	#heatmap_data = gaussian_filter1d(heatmap_data,sigma=sigma_x,axis=1,mode="nearest")

	# Get min/max energies for heatmap.
	if (elim):
		emin_heatmap = elim[0]
		emax_heatmap = elim[1]
	else:
		emin_heatmap = 0
		emax_heatmap = max_energy

	# Heatmap w/ PowerNorm.
	#ax["left"].imshow( heatmap_data,interpolation="bilinear",origin='lower',aspect='auto',cmap='turbo',norm=PowerNorm(gamma=0.12),extent=[0,N_qpts,0,max_energy] )

	# Heatmap w/ LogNorm.
	heatmap_cmap = plt.colormaps["turbo"].copy()
	heatmap_cmap.set_bad("black")
	heatmap_cmap.set_under("black")

	# Heatmap limit scaling.
	vmin = 1E-9 * heatmap_data.max()
	vmax = 1E-3 * heatmap_data.max()
	heatmap_norm = LogNorm(vmin=vmin,vmax=vmax)

	im = ax["left"].imshow( heatmap_data,interpolation="bilinear",origin='lower',aspect='auto',cmap=heatmap_cmap,norm=heatmap_norm,extent=[0,N_qpts,0,max_energy] )

	# Set up colorbar.
	if (qslice_fracs):
		cbar = fig.colorbar(im,ax=ax["left"],location="bottom",fraction=0.04,pad=-0.0025)
	else:
		cbar = fig.colorbar(im,ax=ax["left"],location="bottom",fraction=0.04,pad=0.1000)
	cbar.set_label("vib-EELS Intensity",fontsize=11)

	# Plot HSPs for heatmap.
	if ( len(hsp_locs) > 0 ):
		for i,hsp_loc in enumerate(hsp_locs):
			ax["left"].vlines(hsp_loc,0,max_energy,linestyle="dashed",color="white",alpha=0.5,lw=0.5)
		if (hsp_labels):
			if ( len(hsp_labels) > len(hsp_locs) ):
				print("ERROR: Supplied more HSP labels than there are HSP locations")
				exit()
			ax["left"].set_xticks(hsp_locs[0:len(hsp_labels)],hsp_labels,fontsize=13)
		else:
			hsp_labels = [ x for x in range(len(hsp_locs)) ]
			ax["left"].set_xticks(hsp_locs,hsp_labels,fontsize=13)
	ax["left"].set_xlim(0,hsp_locs[-1])
	ax["left"].set_ylim(emin_heatmap,emax_heatmap)
	ax["left"].yaxis.set_minor_locator(AutoMinorLocator())
	ax["left"].xaxis.set_minor_locator(NullLocator())
	ax["left"].set_ylabel("Energy Loss (meV)",fontsize=12)

	# If we're not doing any q-slices, then just plot this and exit.
	if (not qslice_fracs):
			plt.tight_layout()
			plt.show()
			exit()

	# Plot slice locations on heatmap.
	for i in range(N_slices):
		ax["left"].vlines(slices_todo[i],0,max_energy,color=qslice_colors[i],lw=1.5)

	# Get min/max energies for q-slice plots.
	if (elim):
		emin_qslice = elim[0]
		emax_qslice = elim[1]
	else:
		emin_qslice = 0
		emax_qslice = max_energy * 1.05

	# Now, plot the requested q-slices.
	use_sqrt_J = False
	max_J = -1
	for i in range(N_slices):
		key = f"right_{i}"
		if (use_sqrt_J):
			J_data = np.sqrt(heatmap_data[:,slices_todo[i]])
			ax[key].plot(energies,J_data,color=qslice_colors[i])
			ax[key].set_ylabel(r"$\mathcal{J}^{\frac{1}{2}}$",fontsize=10)
		else:
			J_data = heatmap_data[:,slices_todo[i]]
			ax[key].plot(energies,heatmap_data[:,slices_todo[i]],color=qslice_colors[i])
			ax[key].set_ylabel(r"$\mathcal{J}$",fontsize=10)
		ax[key].set_xlim(emin_qslice,emax_qslice)
		ax[key].ticklabel_format( axis='y',style='sci',scilimits=(0,0) )
		ax[key].minorticks_on()
		# Update maximum J value.
		max_this_J = np.max( J_data.flatten() )
		if (max_this_J > max_J):
			max_J = max_this_J
	# Now ensure that all the y-axes have the same scale, for a meaningful comparison.
	for i in range(N_slices):
		key = f"right_{i}"
		ax[key].set_ylim(0,max_J*1.10)
	plt.xlabel("Energy Loss (meV)",fontsize=12)

	# Finally, show all.
	plt.show()

##########################################################################################
# Low-frequency dielectric tensor.

def plot_eps_lf(seedname, cpts_to_plot, elim=None):

	# New approach to plot layout for eps_lf, with subplot_mosaic.
	# We can hard-code the layouts, since there are only 6 possibilities.
	N_cpts_to_plot = len(cpts_to_plot)
	layouts = { 1: [["a"]],
	2: [["a","a"],
		["b","b"]],
	3: [["a","b"],
		["c","c"]],
	4: [["a","b"],
		["c","d"]],
	5: [["a","b"],
		["c","d"],
		["e","e"]],
	6: [["a","b"],
		["c","d"],
		["e","f"]] }
	sizes = { 1:(6.0,4.5), 2:(6.0,6.5), 3:(8.5,7.0), 4:(8.5,7.0), 5:(8.5,7.0), 6:(8.5,7.0) }
	wspaces = { 1:0.02, 2:0.02, 3:0.04, 4:0.07, 5:0.04, 6:0.08 }
	mosaic = layouts[N_cpts_to_plot]
	fig,axd = plt.subplot_mosaic(mosaic,figsize=sizes[N_cpts_to_plot],constrained_layout=True)
	fig.set_constrained_layout_pads(w_pad=0.05,h_pad=0.05,wspace=wspaces[N_cpts_to_plot],hspace=0.05)
	axs = list(axd.values())

	file = open("./"+seedname+"_eps-lf.dat","r")
	lines = file.readlines()

	units = "inv_cm"

	this_cpt = None
	eps_lf_data = {}
	for line in lines:
		if ("meV" in line):
			units = "meV"
			continue
		if ("Tensor Component" in line):
			if (this_cpt):
				eps_lf_data[this_cpt] = np.array(lfs_real) + 1j*np.array(lfs_imag)
			line = clean_line(line)
			this_cpt = line[-1]
			freqs = []
			lfs_real = []
			lfs_imag = []
			continue
		line = clean_line(line)
		if ( len(line) == 3 ):
			try:
				line = [ float(x) for x in line ]
				freq = line[0]
				lf_real = line[1]
				lf_imag = line[2]
				freqs.append(freq)
				lfs_real.append(lf_real)
				lfs_imag.append(lf_imag)
			except:
				continue

	# Don't forget the last component.
	eps_lf_data[this_cpt] = np.array(lfs_real) + 1j*np.array(lfs_imag)

	for i,cpt in enumerate(cpts_to_plot):
		eps_lf = eps_lf_data[cpt]
		ax = axs[i]
		if (elim):
			emin,emax = elim[0],elim[1]
			ax.set_xlim(emin,emax)
		else:
			ax.set_xlim(0,max(freqs))
		ax.hlines(0,0,max(freqs),color="k",linestyle="dashed",alpha=0.25)
		ax.plot(freqs,eps_lf.real,color="b",label="Real")
		ax.plot(freqs,eps_lf.imag,color="r",label="Imag")
		ax.minorticks_on()
		#ax.set_xticks(fontsize=12)
		#ax.set_yticks(fontsize=12)
		if (units == "meV"):
			ax.set_xlabel(r"Energy (meV)",fontsize=12)
		else:
			ax.set_xlabel(r"Frequency (cm$^{-1}$)",fontsize=12)
		ax.set_ylabel( r"$\epsilon_{_{\mathrm{%s}}}$" % (cpt) , fontsize=16, labelpad=0)
		ax.legend(fontsize=10)

	#plt.tight_layout()
	plt.show()

##########################################################################################
# Polarizability.

def plot_alpha(seedname):

	file = open("./"+seedname+"_alpha.dat","r")
	lines = file.readlines()

	units = "inv_cm"

	all_phis = []
	all_freqs = []
	all_alphas_real = []
	all_alphas_imag = []
	phis = []
	freqs = []
	alphas_real = []
	alphas_imag = []
	current_phi = None
	N_unique_phis = 1
	for line in lines:
		if ("meV" in line):
			units = "meV"
			continue
		if ( ("Phi" in line) or ("Frequency" in line) ):
			continue
		line = clean_line(line)
		if ( len(line) == 4 ):
			try:
				line = [ float(x) for x in line ]
				phi = line[0]
				freq = line[1]
				alpha_real = line[2]
				alpha_imag = line[3]
				if (current_phi is None):
					current_phi = phi
				if ( abs(current_phi - phi) > 1E-6 ):
					all_phis.append(phis)
					all_freqs.append(freqs)
					all_alphas_real.append(alphas_real)
					all_alphas_imag.append(alphas_imag)
					freqs = []
					phis = []
					alphas_real = []
					alphas_imag = []
					current_phi = phi
					N_unique_phis += 1
				phis.append(phi)
				freqs.append(freq)
				alphas_real.append(alpha_real)
				alphas_imag.append(alpha_imag)
			except:
				continue
	# Don't forget last set of values.
	all_phis.append(phis)
	all_freqs.append(freqs)
	all_alphas_real.append(alphas_real)
	all_alphas_imag.append(alphas_imag)
	# Cast to array.
	all_phis = np.array(all_phis)
	all_freqs = np.array(all_freqs)
	all_alphas_real = np.array(all_alphas_real)
	all_alphas_imag = np.array(all_alphas_imag)

	# Just get the unique energies, and unique phis.
	energies = all_freqs[0]
	phis = all_phis[:,0]

	fig = plt.figure(figsize=(13, 5.0), constrained_layout=True)

	gs = fig.add_gridspec(1,6,width_ratios=[2.5, 1.2, 0.08, 2.5, 1.2, 0.08],wspace=0.03)

	ax_re = fig.add_subplot(gs[0, 0], projection="polar")
	ax_im = fig.add_subplot(gs[0, 3], projection="polar")

	gs_re_slices = gs[0, 1].subgridspec(3, 1, hspace=0.25)
	gs_im_slices = gs[0, 4].subgridspec(3, 1, hspace=0.25)

	axes_re = [fig.add_subplot(gs_re_slices[i, 0]) for i in range(3)]
	axes_im = [fig.add_subplot(gs_im_slices[i, 0]) for i in range(3)]

	cax1 = fig.add_subplot(gs[0, 2])
	cax2 = fig.add_subplot(gs[0, 5])

	pcm1 = ax_re.pcolormesh(phis,energies,all_alphas_real.transpose(),shading="auto",cmap="turbo")

	pcm2 = ax_im.pcolormesh(phis,energies,all_alphas_imag.transpose(),shading="auto",cmap="turbo")

	for ax in [ax_re, ax_im]:

		ax.grid(False)

		ax.set_thetamin(-90)
		ax.set_thetamax(90)

		# Angular labels in radians
		ax.set_xticks([-np.pi/2,-np.pi/4,0.0,np.pi/4,np.pi/2])

		ax.set_xticklabels([r"$-\pi/2$",r"$-\pi/4$",r"$0$",r"$\pi/4$",r"$\pi/2$"])

		# Move radial tick labels slightly inward
		ax.set_rlabel_position(-88)

		# Custom radial-axis label
		ax.text(0.12,0.25,"Energy (meV)",rotation=90,va="center",ha="center",transform=ax.transAxes)

	ax_re.set_title(r"$\Re(\alpha)$")
	ax_im.set_title(r"$\Im(\alpha)$")

	slice_phis = [0.0, np.pi/4, np.pi/2]
	slice_labels = [r"$\phi=0$",r"$\phi=\pi/4$",r"$\phi=\pi/2$"]
	slice_colors = ["red","cornflowerblue","orange"]

	rmin = np.nanmin(energies)
	rmax = np.nanmax(energies)

	for phi, label, color, axr, axi in zip(slice_phis, slice_labels, slice_colors, axes_re, axes_im):
		j = np.argmin(np.abs(phis - phi))
		phi_plot = phis[j]

		axr.plot(energies, all_alphas_real[j, :], color=color)
		axi.plot(energies, all_alphas_imag[j, :], color=color)

		axr.set_title(label, fontsize=10, color=color)
		axi.set_title(label, fontsize=10, color=color)

		axr.minorticks_on()
		axi.minorticks_on()

		ax_re.plot([phi_plot, phi_plot], [rmin, rmax],ls="--", lw=1.5, color=color)

		ax_im.plot([phi_plot, phi_plot], [rmin, rmax],ls="--", lw=1.5, color=color)

	for ax in axes_re[:-1] + axes_im[:-1]:
		ax.tick_params(labelbottom=False)

	axes_re[-1].set_xlabel("Energy (meV)")
	axes_im[-1].set_xlabel("Energy (meV)")

	for ax in axes_re:
		ax.set_ylim(np.nanmin(all_alphas_real), np.nanmax(all_alphas_real))

	for ax in axes_im:
		ax.set_ylim(np.nanmin(all_alphas_imag), np.nanmax(all_alphas_imag))

	fig.colorbar(pcm1, cax=cax1)
	fig.colorbar(pcm2, cax=cax2)

	plt.show()

##########################################################################################
# Aloof loss.

def plot_aloof_loss(seedname, elim=None):

	file = open("./"+seedname+"_aloof-loss.dat","r")
	lines = file.readlines()

	units = "inv_cm"

	freqs = []
	losses = []
	for line in lines:
		if ("meV" in line):
			units = "meV"
			continue
		if ( ("Electron" in line) or ("Frequency" in line) or ("Parameter" in line) ):
			continue
		line = clean_line(line)
		if ( len(line) == 2 ):
			try:
				line = [ float(x) for x in line ]
				freq = line[0]
				loss = line[1]
				freqs.append(freq)
				losses.append(loss)
			except:
				continue

	# Get min/max energies for plot.
	if (elim):
		emin,emax = elim[0],elim[1]
	else:
		emin,emax = freqs[0],freqs[-1]

	# Plot aloof loss.
	plt.hlines(0,0,max(freqs),color="k",linestyle="dashed",alpha=0.25)
	plt.plot(freqs,losses,color="b")
	plt.xlim(emin,emax)
	if (units == "inv_cm"):
		plt.xlabel(r"Frequency (cm$^{-1}$)",fontsize=12)
	else:
		plt.xlabel(r"Energy (meV)",fontsize=12)
	plt.ylabel(r"Aloof Loss",fontsize=12)
	plt.tight_layout()
	plt.show()

##########################################################################################
# Oscillator strength tensor.

def plot_ost(seedname):

	file = open("./"+seedname+"_ost.dat","r")
	lines = file.readlines()

	freqs = []
	osc_strength_tensors = []
	N_branches = 0
	for line in lines:
		line = line.replace("\n","")
		line = [ x for x in line.split(" ") if x ]
		if ("frequency:" in line):
			freq = float(line[-2])
			freqs.append(freq)
			N_branches += 1
			continue
		if ( len(line) == 3 ):
			line = [ float(x) for x in line ]
			osc_strength_tensors.append(line)

	# Cast to array, then get trace of OST.
	osc_strength_tensors = np.array(osc_strength_tensors).reshape( (N_branches,3,3) )
	osc_strength_scalars = np.einsum("mii -> m", osc_strength_tensors)

	# Plot.
	plt.vlines(freqs,ymin=0,ymax=osc_strength_scalars,color="b")
	plt.xlabel("Phonon Energy (meV)",fontsize=13)
	plt.ylabel("Oscillator Strength",fontsize=13)
	plt.ylim( 0,max(osc_strength_scalars)*1.10 )
	plt.xlim( 0,max(freqs)*1.10 )
	plt.xticks(fontsize=11)
	plt.yticks(fontsize=11)
	plt.minorticks_on()

	plt.show()
	plt.tight_layout()

##########################################################################################
# Aloof w/ OST.

def plot_aloof_ost(seedname, elim=None):

	file = open("./"+seedname+"_aloof-loss.dat","r")
	lines = file.readlines()

	units = "inv_cm"

	aloof_freqs = []
	losses = []
	for line in lines:
		if ("meV" in line):
			units = "meV"
			continue
		if ( ("Electron" in line) or ("Frequency" in line) ):
			continue
		line = clean_line(line)
		if ( len(line) == 2 ):
			try:
				line = [ float(x) for x in line ]
				freq = line[0]
				loss = line[1]
				aloof_freqs.append(freq)
				losses.append(loss)
			except:
				continue

	fig,axs = plt.subplots(2,1,sharex=True)

	# Get min/max energies for plot.
	if (elim):
		emin,emax = elim[0],elim[1]
	else:
		emin,emax = 0,max(aloof_freqs)

	# Plot aloof loss.
	axs[0].hlines(0,0,max(aloof_freqs),color="k",linestyle="dashed",alpha=0.25)
	axs[0].plot(aloof_freqs,losses,color="k",lw=0.85)
	axs[0].set_xlim(emin,emax)
	axs[0].set_ylim( 0,max(losses)*1.10 )
	axs[0].set_ylabel("Aloof Loss",fontsize=12)
	plt.yticks(fontsize=11)
	axs[0].minorticks_on()

	# Now, plot oscillator strength tensor.
	file = open("./"+seedname+"_ost.dat","r")
	lines = file.readlines()

	ost_freqs = []
	osc_strength_tensors = []
	N_branches = 0
	for line in lines:
			line = line.replace("\n","")
			line = [ x for x in line.split(" ") if x ]
			if ("frequency:" in line):
					freq = float(line[-2])
					ost_freqs.append(freq)
					N_branches += 1
					continue
			if ( len(line) == 3 ):
					line = [ float(x) for x in line ]
					osc_strength_tensors.append(line)

	# Cast to array, then get trace of OST.
	osc_strength_tensors = np.array(osc_strength_tensors).reshape( (N_branches,3,3) )
	osc_strength_scalars = np.einsum("mii -> m", osc_strength_tensors)

	# Normalize OST.
	osc_strength_scalars = osc_strength_scalars / max( osc_strength_scalars.flatten() )

	# Plot.
	axs[1].vlines(ost_freqs,ymin=0,ymax=osc_strength_scalars,color="b")
	axs[1].set_xlabel("Energy Loss (meV)",fontsize=12)
	axs[1].set_ylabel(r"$\mathrm{Tr}[S^{m}_{\alpha,\beta}]$",fontsize=12)
	axs[1].set_ylim( 0,max(osc_strength_scalars)*1.10 )
	plt.xticks(fontsize=11)
	plt.yticks(fontsize=11)
	axs[1].minorticks_on()

	fig.subplots_adjust(hspace=0)

	plt.show()

##########################################################################################
# Atomic form factor.

def plot_aff(seedname, hsp_labels=None):

	# Plot an aff curve along the high symmetry line, for each element.
	aff_max = 0
	aff_min = 0
	file = open("./"+seedname+"_aff.dat","r")
	lines = file.readlines()

	element_types = []
	affs = []
	read_hsp_locs = False
	read_data = False
	for line in lines:
		if ("High-symmetry points" in line):
			read_hsp_locs = True
			continue
		if (read_hsp_locs):
			line = clean_line(line)
			hsp_locs = [ int(x) - 1 for x in line ]
			read_hsp_locs = False
			continue
		if ("aff(|q|)" in line):
			read_data = True
			continue
		if (read_data):
			try:
				line = line.replace("\n","")
				line = [ x for x in line.split(" ") if x ]
				if ( len(line) == 0 ):
					continue
				el_type = str(line[0])
				aff = float(line[-1])
				affs.append(aff)
				if (el_type not in element_types):
					element_types.append(el_type)
			except:
				continue
	N_elements = len(element_types)
	N_qpts = int( len(affs)/N_elements )
	affs = np.array(affs).reshape( (N_elements,N_qpts) )
	for i,el_string in enumerate(element_types):
		plt.plot(affs[i,:],label=el_string)
	aff_min = min(affs.flatten())
	aff_max = max(affs.flatten())

	# HSPs.
	if ( len(hsp_locs) > 0 ):
		for i,hsp_loc in enumerate(hsp_locs):
			plt.vlines(hsp_loc,0,aff_max*1.05,linestyle="dashed",color="k",alpha=0.5,lw=0.5)
		if (hsp_labels):
			if ( len(hsp_labels) > len(hsp_locs) ):
				print("ERROR: Supplied more HSP labels than there are HSP locations")
				exit()
			plt.xticks(hsp_locs[0:len(hsp_labels)],hsp_labels,fontsize=13)
		else:
			hsp_labels = [ x for x in range(len(hsp_locs)) ]
			plt.xticks(hsp_locs,hsp_labels,fontsize=13)

	plt.xlim(0,hsp_locs[-1])
	plt.ylim(aff_min*0.90,aff_max*1.10)
	plt.yticks(fontsize=12)
	plt.ylabel(r"$\mathrm{AFF}(\mathrm{\textbf{q}})$",fontsize=12)
	plt.legend(ncols=2,fontsize=10)
	plt.tight_layout()
	plt.show()

##########################################################################################
# Debye-Waller factor.

def plot_dwf(seedname, hsp_labels=None):

	# Plot a DWF curve along the high symmetry line, for each element.
	dwf_max = 1
	dwf_min = 1

	file = open("./"+seedname+"_dwf.dat","r")
	lines = file.readlines()

	atom_labels = []
	dwfs = []
	read_hsp_locs = False
	read_data = False
	for line in lines:
		if ("High-symmetry points" in line):
			read_hsp_locs = True
			continue
		if (read_hsp_locs):
			line = clean_line(line)
			hsp_locs = [ int(x) - 1 for x in line ]
			read_hsp_locs = False
			continue
		if ("dwf(q)" in line):
			read_data = True
			continue
		if (read_data):
			try:
				line = line.replace("\n","")
				line = [ x for x in line.split(" ") if x ]
				if ( len(line) == 0 ):
					continue
				atom_label = str(line[0])
				dwf = float(line[-1])
				dwfs.append(dwf)
				if (atom_label not in atom_labels):
					atom_labels.append(atom_label)
			except:
				continue
	N_atoms = len(atom_labels)
	N_qpts = int( len(dwfs) / N_atoms )
	dwfs = np.array(dwfs).reshape( (N_atoms,N_qpts) )
	for i,label_string in enumerate(atom_labels):
		plt.plot(dwfs[i,:],label=label_string)
	dwf_min = min(dwfs.flatten())
	dwf_max = max(dwfs.flatten())

	# Plot horizontal reference line at DWF=1
	plt.hlines(1,0,hsp_locs[-1],linestyle="dashed",alpha=0.20,color="k")

	# HSPs.
	if ( len(hsp_locs) > 0 ):
		for i,hsp_loc in enumerate(hsp_locs):
			plt.vlines(hsp_loc,0,dwf_max*1.05,linestyle="dashed",color="k",alpha=0.5,lw=0.5)
		if (hsp_labels):
			if ( len(hsp_labels) > len(hsp_locs) ):
				print("ERROR: Supplied more HSP labels than there are HSP locations")
				exit()
			plt.xticks(hsp_locs[0:len(hsp_labels)],hsp_labels,fontsize=13)
		else:
			hsp_labels = [ x for x in range(len(hsp_locs)) ]
			plt.xticks(hsp_locs,hsp_labels,fontsize=13)

	plt.xlim(0,hsp_locs[-1])
	plt.ylim(dwf_min*0.995,dwf_max*1.005)
	plt.yticks(fontsize=11)
	plt.ylabel(r"$\mathrm{DWF}(\textbf{q})$",fontsize=12)
	plt.legend(ncols=2,fontsize=10)
	plt.tight_layout()
	plt.show()

##########################################################################################
# Impact-EELS intensity.
#   --> This uses a cmap for the impact EELS-intensity, with a skeleton dispersion
#       plotted in light grey underneath. This is excellent for displaying the
#       intensity information, but has the disadvantage that some band-crossings may
#       be harder to see.

def plot_impact_intensity(seedname, hsp_labels=None, elim=None):

	file = open("./"+seedname+"_impact-intensity.dat","r")
	lines = file.readlines()

	units = "inv_cm"

	N_bnd = 0
	freqs = []
	eels_intensities = []
	read_hsp_locs = False
	read_data = False
	for line in lines:
		if ("meV" in line):
			units = "meV"
			continue
		if ("Evaluated" in line):
			line = clean_line(line)
			qa_ref,qb_ref,qc_ref = float(line[-3]),float(line[-2]),float(line[-1])
			continue
		if ("High-symmetry points" in line):
			read_hsp_locs = True
			continue
		if ("q_idx" in line):
			continue
		if (read_hsp_locs):
			line = clean_line(line)
			hsp_locs = [ int(x) - 1 for x in line ]
			read_hsp_locs = False
			read_data = True
			continue
		if (read_data):
			try:
				line = clean_line(line)
				if ( len(line) == 0 ):
					continue
				qpt_idx = int( line[0] )
				bnd_idx = int( line[1] )
				freq = float( line[2] )
				eels_intensity = float( line[3] )
				freqs.append(freq)
				eels_intensities.append(eels_intensity)
				if (bnd_idx > N_bnd):
					N_bnd = bnd_idx
			except:
				continue

	# Last qpt_idx gives total number of q-points.
	N_qpts = qpt_idx

	# Cast to array.
	freqs = np.array(freqs).reshape( (N_qpts,N_bnd) )
	eels_intensities = np.array(eels_intensities).reshape( (N_qpts,N_bnd) )

	# Get maximum energy.
	max_energy = np.max( freqs.flatten() )

	# Scaling factor for EELS intensity scatter points.
	SCTR_SCALE = 75 / np.max(eels_intensities.flatten()) # 50

	# Normalize colour scale between lower and upper percentiles.
	I_all = eels_intensities.flatten()
	#perc_min = 20.0
	#perc_max = 99.5
	perc_min = 30.0  # 23.0
	perc_max = 99.85 # 99.7
	vmin = np.percentile(I_all,perc_min)
	vmax = np.percentile(I_all,perc_max)

	# Sometimes, at the lowest percentile, vmin can be zero, which messes up the LogNorm.
	# So slowly bump it up until we hit a non-zero vmin.
	i = 1
	while (vmin < 1E-24):
		vmin = np.percentile(I_all, perc_min + 3*i)
		i += 1
	norm = LogNorm( vmin=vmin,vmax=vmax )

	# Helper function for EELS-intensity plot.
	def colored_line(k, omega, intensity, cmap, norm, lw=2):
		points = np.array([k, omega]).T.reshape(-1, 1, 2)
		segments = np.concatenate([points[:-1], points[1:]], axis=1)
		I_seg = 0.5 * (intensity[:-1] + intensity[1:])
		lc = LineCollection(segments, cmap=cmap, norm=norm)
		lc.set_array(I_seg)
		lc.set_linewidth(lw)
		return lc

	# Helper function to truncate cmap range.
	def truncate_colormap(cmap, minval=0.05, maxval=0.95, n=256):
		#return mcolors.LinearSegmentedColormap.from_list(f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',cmap(np.linspace(minval, maxval, n)))
		return LinearSegmentedColormap.from_list(f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',cmap(np.linspace(minval, maxval, n)))

	# Set up figure.
	fig,ax = plt.subplots()
	cmap = plt.get_cmap("turbo")
	cmap = truncate_colormap(cmap, 0.03, 0.88)

	# Loop over bands.
	# Plot skeleton for dispersion in grey.
	# EELS-active bands are plotted on a varying colour scale.
	for bnd in range(N_bnd):
		qpts = [ x for x in range(N_qpts) ]
		ax.plot(qpts,freqs[:,bnd],color="lightgray",lw=0.8,zorder=1)
		lc = colored_line(qpts,freqs[:,bnd],eels_intensities[:,bnd],cmap=cmap,norm=norm,lw=2)
		ax.add_collection(lc)

	# Colorbar.
	sm = cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])
	cbar = plt.colorbar(sm, ax=ax)
	cbar.set_label("vib-EELS Intensity",fontsize=11)

	# Get min/max energies for plot.
	if (elim):
		emin_plot = elim[0]
		emax_plot = elim[1]
	else:
		emin_plot = 0
		emax_plot = max_energy * 1.05

	# HSPs.
	if ( len(hsp_locs) > 0 ):
		for i,hsp_loc in enumerate(hsp_locs):
			plt.vlines(hsp_loc,emin_plot,emax_plot,linestyle="dashed",color="k",alpha=0.5,lw=0.5)
		if (hsp_labels):
			if ( len(hsp_labels) > len(hsp_locs) ):
				print("ERROR: Supplied more HSP labels than there are HSP locations")
				exit()
			plt.xticks(hsp_locs[0:len(hsp_labels)],hsp_labels,fontsize=13)
		else:
			hsp_labels = [ x for x in range(len(hsp_locs)) ]
			plt.xticks(hsp_locs,hsp_labels,fontsize=13)

	plt.xlim(0,hsp_locs[-1])
	if (elim):
		plt.ylim(elim[0],elim[1])
	else:
		plt.ylim(emin_plot,emax_plot)
	if (units == "meV"):
		plt.ylabel(r"Phonon Energy (meV)",fontsize=12)
	else:
		plt.ylabel(r"Phonon Frequency (cm$^{-1}$)",fontsize=12)
	plt.gca().xaxis.set_minor_locator(NullLocator())
	plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
	plt.tight_layout()
	plt.show()

##########################################################################################
# Just plot phonon dispersion.

def plot_dispersion(seedname, hsp_labels=None, elim=None):

	file = open("./"+seedname+"_impact-intensity.dat","r")
	lines = file.readlines()

	units = "inv_cm"

	N_bnd = 0
	freqs = []
	read_hsp_locs = False
	read_data = False
	for line in lines:
		if ("meV" in line):
			units = "meV"
			continue
		if ("Evaluated" in line):
			line = clean_line(line)
			qa_ref,qb_ref,qc_ref = float(line[-3]),float(line[-2]),float(line[-2])
			continue
		if ("High-symmetry points" in line):
			read_hsp_locs = True
			continue
		if ("q_idx" in line):
			continue
		if (read_hsp_locs):
			line = clean_line(line)
			hsp_locs = [ int(x) - 1 for x in line ]
			read_hsp_locs = False
			read_data = True
			continue
		if (read_data):
			try:
				line = clean_line(line)
				if ( len(line) == 0 ):
					continue
				qpt_idx = int( line[0] )
				bnd_idx = int( line[1] )
				freq = float( line[2] )
				freqs.append(freq)
				if (bnd_idx > N_bnd):
					N_bnd = bnd_idx
			except:
				continue

	# Last qpt_idx gives total number of q-points.
	N_qpts = qpt_idx

	# Cast to array.
	freqs = np.array(freqs).reshape( (N_qpts,N_bnd) )

	# Get maximum energy.
	max_energy = np.max( freqs.flatten() )

	for bnd in range(N_bnd):
		qpts = [ x for x in range(N_qpts) ]

		# Plot each band with a different colour.
		plt.plot(freqs[:,bnd],lw=1,zorder=1,rasterized=True)

	# Get min/max energies for plot.
	if (elim):
		emin_plot = elim[0]
		emax_plot = elim[1]
	else:
		emin_plot = 0
		emax_plot = max_energy * 1.05

	# HSPs.
	if ( len(hsp_locs) > 0 ):
		for i,hsp_loc in enumerate(hsp_locs):
			plt.vlines(hsp_loc,emin_plot,emax_plot,linestyle="dashed",color="k",alpha=0.5,lw=0.5)
		if (hsp_labels):
			if ( len(hsp_labels) > len(hsp_locs) ):
				print("ERROR: Supplied more HSP labels than there are HSP locations")
				exit()
			plt.xticks(hsp_locs[0:len(hsp_labels)],hsp_labels,fontsize=13)
		else:
			hsp_labels = [ x for x in range(len(hsp_locs)) ]
			plt.xticks(hsp_locs,hsp_labels,fontsize=13)

	plt.xlim(0,hsp_locs[-1])
	plt.ylim(emin_plot,emax_plot)
	if (units == "meV"):
		plt.ylabel(r"Phonon Energy (meV)",fontsize=12)
	else:
		plt.ylabel(r"Phonon Frequency (cm$^{-1}$)",fontsize=12)
	plt.gca().xaxis.set_minor_locator(NullLocator())
	plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
	plt.tight_layout()
	plt.show()

##########################################################################################
# Run.

parser = argparse.ArgumentParser()
parser.add_argument("plot_todo", type=str, help="Type of plot to produce. Should be one of: eps_lf, alpha, aff, dwf, "
"impact_intensity, dispersion, aloof_loss, osc_strength, aloof_ost, or impact_heatmap")

# HSPs.
parser.add_argument("-hsp_labels", type=str, nargs="+", help="High-symmetry labels, as a space separated list (optional)")
parser.add_argument("-qslices", type=float, nargs="+", help="Points along the path (as a fraction of total path "
"length) at which to draw vertical slices; optional, and only relevant for the 'impact_heatmap' plot type")
parser.add_argument("-cpts", type=str, nargs="+", help="Tensor components to plot for eps_lf (default is just 'xx')")
parser.add_argument("-elim", type=float, nargs=2, help="Minimum / Maximum energies for plot")

# Parse arguments.
args = parser.parse_args()
plot_todo = args.plot_todo
hsp_labels = args.hsp_labels
qslices = args.qslices
eps_lf_cpts = args.cpts
elim = args.elim

# Render special points appropriately w/ LaTeX.
if (hsp_labels):
	for i,lbl in enumerate(hsp_labels):
		# All 'G' labels are Gamma-points.
		if (lbl == "G"):
			hsp_labels[i] = r"$\Gamma$"
		# The presence of 'b' in a label (e.g. AbX) denotes a path break.
		if ("b" in lbl):
			hsp_labels[i] = r"%s$\vert$%s" % (lbl[0],lbl[-1])

# Default for eps_lf tensor component, if nothing was specified.
if (not eps_lf_cpts):
	eps_lf_cpts = ["xx"]

##########################################################################################
# Execute relevant plot.

# Get seedname from .eels_data file.
eels_file = glob.glob("*.eels_data")
if (not eels_file):
	eels_file = glob.glob("*.eels_input")
if (not eels_file):
	eels_file = glob.glob("*.odd")
seedname = eels_file[0].split(".")[0]

if (plot_todo == "eps_lf"):
	plot_eps_lf(seedname, cpts_to_plot=eps_lf_cpts, elim=elim)
elif (plot_todo == "alpha"):
	plot_alpha(seedname)
elif (plot_todo == "aff"):
	plot_aff(seedname, hsp_labels=hsp_labels)
elif (plot_todo == "dwf"):
	plot_dwf(seedname, hsp_labels=hsp_labels)
elif (plot_todo == "impact_intensity"):
	plot_impact_intensity(seedname, hsp_labels=hsp_labels, elim=elim)
elif (plot_todo == "dispersion"):
	plot_dispersion(seedname, hsp_labels=hsp_labels, elim=elim)
elif (plot_todo == "aloof_loss"):
	plot_aloof_loss(seedname, elim=elim)
elif (plot_todo == "ost"):
	plot_ost(seedname)
elif (plot_todo == "aloof_ost"):
	plot_aloof_ost(seedname, elim=elim)
elif (plot_todo == "impact_heatmap"):
	plot_impact_heatmap(seedname, hsp_labels=hsp_labels, qslice_fracs=qslices, elim=elim)
else:
	print("ERROR: Plot type not recognised")
	print("Should be one of: eps_lf, alpha, aff, dwf, impact_intensity, dispersion, aloof_loss, "
	"ost, aloof_ost, or impact_heatmap")
	exit()
