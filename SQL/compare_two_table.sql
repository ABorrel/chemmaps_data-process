SELECT id,name
FROM(

SELECT public.chem_descriptor_opera_name_new.id , public.chem_descriptor_opera_name_new.name
FROM public.chem_descriptor_opera_name_new
UNION ALL
SELECT update_public.chem_descriptor_opera_name.id , update_public.chem_descriptor_opera_name.name
FROM update_public.chem_descriptor_opera_name) tbl
GROUP BY id, name
HAVING count(*) = 1
ORDER BY id;