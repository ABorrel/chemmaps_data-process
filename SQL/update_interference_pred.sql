-- need to be in the schema to update

update update_public.chemical_description 
set desc_opera = public.chemical_description.desc_opera 
from public.chemical_description
where update_public.chemical_description.inchikey = public.chemical_description.inchikey 