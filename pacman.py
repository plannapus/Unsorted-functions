# PacMan Trimming (Lazarus et al. 2012)
# neptune_object is a list of dictionaries (typically the result of a postgreSQL query using psycopg2)
# Each dictionary should at least contains a taxon_id field (giving a taxon) and a sample_age_ma (giving an age) field.
# pactop and pacbottom are respectively the top and bottom percentage of the range to trim.

def PacMan(neptune_object, pactop, pacbot):
    from math import floor
    id = 'taxon_id'
    age_id = 'sample_age_ma'   
    species = set([k[id] for k in neptune_object])
    res = []
    for i in species:
        species_occurrences = [k for k in neptune_object if k[id] == i]
        species_occurrences = sorted(species_occurrences, key=lambda k:k[age_id])
        Nocc = len(species_occurrences)
        nb_top = int(floor(Nocc * pactop /100))
        nb_bot = int(floor(Nocc * pacbot /100))
        trimmed = species_occurrences[nb_top:(Nocc-nb_bot-1)]
        for j in trimmed: res.append(j)
    
    res = sorted(res, key=lambda k:k[age_id])
    return(res)
