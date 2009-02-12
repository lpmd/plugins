#
#
#

EXTRA_DIST = update-svnrev.sh

version.h:
	./update-svnrev.sh .

all:
	$(MAKE) version.h
	$(MAKE) all-am

install:
	./update-svnrev.sh .
	$(MAKE) install-am

install-exec:
	./update-svnrev.sh .
	$(MAKE) install-exec-am

CLEANFILES = version.h

noinst_LTLIBRARIES = plugincommon.la
plugincommon_la_SOURCES = plugincommon.h plugincommon.cc version.h

#
#
#
plugindir = $(libdir)/lpmd

AM_CPPFLAGS = -Wall -pedantic
AM_LDFLAGS = -module -avoid-version plugincommon.lo

#
# plugin_LTLIBRARIES contiene la lista de plugins de liblpmd (como librerias de libtool)
#
# plugin_LTLIBRARIES contains the list of liblpmd plugins (as libtool libraries)
#
plugin_LTLIBRARIES =    angdist.la beeman.la berendsen.la buckingham.la cell.la \
                        cellscaling.la cna.la constantforce.la cordnum.la cordnumfunc.la \
                        crystal2d.la crystalbcc.la crystalfcc.la crystalhcp.la crystalsc.la \
                        densityprofile.la displace.la dispvol.la dlpoly.la energy.la \
                        euler.la ewald.la extravel.la fastlj.la finnissinclair.la gdr.la \
                        gupta.la harmonic.la leapfrog.la lennardjones.la linkedcell.la lpmd.la \
                        localpressure.la minimumimage.la mol2.la moleculecm.la morse.la msd.la \
                        nosehoover.la nullintegrator.la nullpairpotential.la nullmetalpotential.la \
                        nullpotential.la pairdistances.la povray.la pressure.la printatoms.la \
                        quenchedmd.la rawbinary.la replicate.la rotate.la rvcorr.la selectatoms.la shear.la \
                        skewstart.la suttonchen.la thermalneedle.la temperature.la tempprofile.la \
                        tempscaling.la vacf.la vasp.la veldist.la velocityverlet.la verlet.la \
                        verletlist.la voronoi.la xyz.la zlp.la

#
#
#

zlp_la_LDFLAGS = -module -avoid-version -lz

#
#
#

angdist_la_SOURCES = angdist.h angdist.cc 
beeman_la_SOURCES = beeman.h beeman.cc
berendsen_la_SOURCES = berendsen.h berendsen.cc
buckingham_la_SOURCES = buckingham.h buckingham.cc
cell_la_SOURCES = cell.h cell.cc
cellscaling_la_SOURCES = cellscaling.h cellscaling.cc
cna_la_SOURCES = cna.h cna.cc
constantforce_la_SOURCES = constantforce.h constantforce.cc
cordnum_la_SOURCES = cordnum.h cordnum.cc
cordnumfunc_la_SOURCES = cordnumfunc.h cordnumfunc.cc
crystal2d_la_SOURCES = crystal2d.h crystal2d.cc
crystalbcc_la_SOURCES = crystalbcc.h crystalbcc.cc
crystalfcc_la_SOURCES = crystalfcc.h crystalfcc.cc
crystalhcp_la_SOURCES = crystalhcp.h crystalhcp.cc
crystalsc_la_SOURCES = crystalsc.h crystalsc.cc
densityprofile_la_SOURCES = densityprofile.h densityprofile.cc
displace_la_SOURCES = displace.h displace.cc
dispvol_la_SOURCES = dispvol.h dispvol.cc
dlpoly_la_SOURCES = dlpoly.h dlpoly.cc
energy_la_SOURCES = energy.h energy.cc
euler_la_SOURCES = euler.h euler.cc
ewald_la_SOURCES = ewald.h ewald.cc
extravel_la_SOURCES = extravel.h extravel.cc
fastlj_la_SOURCES = fastlj.h fastlj.cc
finnissinclair_la_SOURCES = finnissinclair.h finnissinclair.cc
gdr_la_SOURCES = gdr.h gdr.cc
gupta_la_SOURCES = gupta.h gupta.cc
harmonic_la_SOURCES = harmonic.h harmonic.cc
leapfrog_la_SOURCES = leapfrog.h leapfrog.cc
lennardjones_la_SOURCES = lennardjones.h lennardjones.cc
linkedcell_la_SOURCES = linkedcell.h linkedcell.cc
localpressure_la_SOURCES = localpressure.h localpressure.cc
lpmd_la_SOURCES = lpmd.h lpmd.cc
minimumimage_la_SOURCES = minimumimage.h minimumimage.cc
mol2_la_SOURCES = mol2.h mol2.cc
moleculecm_la_SOURCES = moleculecm.h moleculecm.cc
morse_la_SOURCES = morse.h morse.cc
msd_la_SOURCES = msd.h msd.cc
nosehoover_la_SOURCES = nosehoover.h nosehoover.cc
nullintegrator_la_SOURCES = nullintegrator.h nullintegrator.cc
nullpairpotential_la_SOURCES = nullpairpotential.h nullpairpotential.cc
nullmetalpotential_la_SOURCES = nullmetalpotential.h nullmetalpotential.cc
nullpotential_la_SOURCES = nullpotential.h nullpotential.cc
pairdistances_la_SOURCES = pairdistances.h pairdistances.cc
povray_la_SOURCES = povray.h povray.cc
pressure_la_SOURCES = pressure.h pressure.cc
printatoms_la_SOURCES = printatoms.h printatoms.cc
quenchedmd_la_SOURCES = quenchedmd.h quenchedmd.cc
rawbinary_la_SOURCES = rawbinary.h rawbinary.cc
replicate_la_SOURCES = replicate.h replicate.cc
rotate_la_SOURCES = rotate.h rotate.cc
rvcorr_la_SOURCES = rvcorr.h rvcorr.cc
selectatoms_la_SOURCES = selectatoms.h selectatoms.cc
shear_la_SOURCES = shear.h shear.cc
skewstart_la_SOURCES = skewstart.h skewstart.cc
suttonchen_la_SOURCES = suttonchen.h suttonchen.cc
thermalneedle_la_SOURCES = thermalneedle.h thermalneedle.cc
temperature_la_SOURCES = temperature.h temperature.cc
tempprofile_la_SOURCES = tempprofile.h tempprofile.cc
tempscaling_la_SOURCES = tempscaling.h tempscaling.cc
vacf_la_SOURCES = vacf.h vacf.cc
vasp_la_SOURCES = vasp.h vasp.cc
veldist_la_SOURCES = veldist.h veldist.cc
velocityverlet_la_SOURCES = velocityverlet.h velocityverlet.cc
verlet_la_SOURCES = verlet.h verlet.cc
verletlist_la_SOURCES = verletlist.h verletlist.cc
voronoi_la_SOURCES = voronoi.h voronoi.cc
xyz_la_SOURCES = xyz.h xyz.cc
zlp_la_SOURCES = zlp.h zlp.cc

