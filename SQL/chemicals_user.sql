-- Table: public.chemicals_user

-- DROP TABLE public.chemicals_user;

CREATE TABLE public.chemicals_user
(
    id bigint,
    smiles_origin text COLLATE pg_catalog."default",
    smiles_clean text COLLATE pg_catalog."default",
    inchikey text COLLATE pg_catalog."default",
    dsstox_id text COLLATE pg_catalog."default",
    drugbank_id text COLLATE pg_catalog."default",
    mol_clean mol,
    casn text COLLATE pg_catalog."default",
    name text COLLATE pg_catalog."default",
    status text COLLATE pg_catalog."default"
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chemicals_user
    OWNER to postgres;