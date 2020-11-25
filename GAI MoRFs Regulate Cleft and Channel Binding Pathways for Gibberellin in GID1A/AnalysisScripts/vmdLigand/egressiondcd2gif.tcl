# name of an input dcd file
set dcd_name  [lindex $argv 0]
set out_name  [lindex $argv 1]
set freq      [lindex $argv 2]
set name $dcd_name
mol load  parm7 ./ref.prmtop dcd $name
set mol1 [molinfo top get id]
put "END INI"
set outfile $out_name-movie.gif
#put "$outfile"

package require pbctools
pbc set [pbc get -all] -all
pbc wrap -center com -centersel protein -compound res -all

set nf [molinfo $mol1 get numframes]

scale by 18
color Display Background white
#mol modselect 1 top protein
mol modselect 0 top residue 1 to 110
mol modstyle 0 top NewCartoon
mol modcolor 0 top name
mol addrep top                                                                                                        
mol modselect 1 top {resname "UNL"}
mol modcolor 1 top timestep                
mol modstyle 1 top vdw
mol addrep top
mol modselect 2 top residue 111 to 440
mol modstyle 2 top NewCartoon
mol modcolor 2 top colorid 3
rotate x by 180
rotate y by 90
for {set i 0} {$i < $nf} {incr i $freq} {
	set filename snap.[format "%04d" $i].tga
	animate goto $i
	display update
	render Tachyon $filename "/apps/pkg/vmd/1.9.3/text/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 -trans_vmd -mediumshade %s -res 1024 1024 -format TARGA -o %s
	rotate y by 1
}
exec convert -delay 10 -depth 8 -loop 0 tga:snap.*.tga $outfile
exec rm snap*
quit

#proc make_trajectory_movie {freq} {
#	# get the number of frames in the movie
#	set num [molinfo top get numframes]
#	# loop through the frames
#	for {set i 0} {$i < $num} {incr i $freq} {
#		# go to the given frame
#		animate goto $i
#               # for the display to update
#                display update
#		# take the picture
#		set filename snap.[format "%04d" [expr $i/$freq]].rgb
#		render snapshot $filename
#	}
#}
