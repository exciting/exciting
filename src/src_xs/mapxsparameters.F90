
! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mapxsparameters
	use modmain
	use modxs
	implicit none
	nosym=nosymxs
	ngridk(:)=ngridkxs(:)
	reducek=reducekxs
	vkloff(:)=vkloffxs(:)
	ngridq(:)=ngridqxs(:)
	reduceq=reduceqxs
	rgkmax=rgkmaxxs
	swidth=swidthxs
	lmaxapw=lmaxapwxs
	lmaxmat=lmaxmatxs
	nempty=nemptyxs
end subroutine
