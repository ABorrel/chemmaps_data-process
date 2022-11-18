-- View: public.mvwchemmap_chemcials

-- DROP MATERIALIZED VIEW public.mvwchemmap_chemcials;

CREATE MATERIALIZED VIEW public.mvwchemmap_chemcials
TABLESPACE pg_default
AS
 SELECT chem.smiles_origin,
    chem.smiles_clean,
    chem.inchikey,
    chem.dsstox_id,
    chem.drugbank_id,
    info.tox21,
    info.pfas
   FROM chemicals chem
     LEFT JOIN chemical_info info ON chem.id = info.chem_id
WITH DATA;

ALTER TABLE public.mvwchemmap_chemcials
    OWNER TO postgres;