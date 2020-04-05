a=$1
for i in serine threonine cysteine valine leucine isoleucine methionine phenylalanine tyrosine tryptophan aspartic glutamic asparagine glutamine histidine_basic histidine_neutral lysine arginine nma 

do	
	python KB_integral_processing.py KB_integral_analysis/rdf-cooranal-urea_"$i"_"$a".dat KB_integral_analysis/rdf-cooranal-wate_"$i"_"$a".dat "$i"_"$a"murea.pdb "$i"
	#python variance_energy.py "$i"_te.dat "$i"_ele.dat "$i"_vdw.dat "$i"

done
