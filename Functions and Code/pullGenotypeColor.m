
function kolor = pullGenotypeColor(geno)

switch geno
    case 'Berlin-WT'
        c = 'Purple';
    case 'CantonS'
        c = 'Teal';
    case 'Poxn_mutant'
        c = 'Red';
    case 'Orco_mutant'
        c = 'Yellow';
    case 'IsoD1'
        c = 'Dodgerblue';
    case 'W(1118)'
        c = 'Moccasin';

end
       
kolor = Color(c);