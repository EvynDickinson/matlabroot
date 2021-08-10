

rental_costs = (1600:100:2200);
rental_duration = (9:12);
august_rent = 850*2;
sub_august_rent = 675*2;
houserent = 2200;
move_costs = 700;

% condition if we stayed in the house and Sarah and Conor moved out early
stay.august = 675*2;
stay.sep_june = stay.august + (9*houserent);
stay.sep_july = stay.august + (10*houserent);
stay.sep_aug = stay.august + (11*houserent);
stay.sep_sep = stay.august + (12*houserent);

for irent = 1:length(rental_costs)
    rent = rental_costs(irent);
    % condition if we started renting a new place in august
    augstart(irent).august = 850*2 + rent + move_costs;
    augstart(irent).aug_june = augstart(irent).august + (10*rent);
    augstart(irent).aug_july = augstart(irent).august + (11*rent);
    augstart(irent).aug_aug = augstart(irent).august + (12*rent);

    % condition if we started renting a new place in august
    sepstart(irent).august = 850*2 + move_costs;
    sepstart(irent).sep_june = sepstart(irent).august + (9*rent);
    sepstart(irent).sep_july = sepstart(irent).august + (10*rent);
    sepstart(irent).sep_aug = sepstart(irent).august + (11*rent);
    sepstart(irent).sep_sep = sepstart(irent).august + (12*rent);
end

buff = 1500;
stay.sep_june + buff
stay.sep_july + buff
stay.sep_aug + buff
stay.sep_sep + buff

stay.sep_june - buff
stay.sep_july - buff
stay.sep_aug - buff
stay.sep_sep - buff



