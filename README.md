    # UjiLity - UJI muLtiple prImer uTilitY
    A script to design sets of oligonucleotide primers with controlled
    degeneracy aimed at amplifying marker genes with very high sequence
    diversity (e.g. virus diversity studies).
    
    The script was written in the spring of 2015 in Uji, Japan, by 
    Pascal Hingamp as a Kyoto University visiting scholar in the lab
    of Hiro Ogata whose team members provided instrumental insights during
    development: Hiro Ogata, Takashi Yoshida, Susumu Goto, Tomoko Mihara
    & Yosuke Nishimura. The resulting primer design strategy was validated 
    by an experimental test of a UjiLity designed set of primers (called
    MEGAPRIMER) that demonstrated detection of giant viruses diversity
    in marine water: https://www.mdpi.com/1999-4915/10/9/496
    
    ## Version 2.17 (2022-05-28) pascal.hingamp@univ-amu.fr
    		
    Warning: as often, this script was never meant to grow so large.
    It therefore rather monolithic, severely lacking in structure,
    as well as begging for richer comments... It does however do the job.
    The STDOUT output is very terse, basically showing progress, whilst
    the STDERR output is over verbose and should be redirected to a file.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
