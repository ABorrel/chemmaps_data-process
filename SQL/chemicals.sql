-- Table: public.chemicals

-- DROP TABLE chemicals;
create sequence chem_seq;

CREATE TABLE chemicals
(
    id bigint NOT NULL DEFAULT nextval('chem_seq'::regclass),
    smiles_origin text COLLATE pg_catalog."default",
    smiles_clean text COLLATE pg_catalog."default",
    inchikey text COLLATE pg_catalog."default",
    dsstox_id text COLLATE pg_catalog."default",
    drugbank_id text COLLATE pg_catalog."default",
    casn text COLLATE pg_catalog."default",
    name text COLLATE pg_catalog."default",
    CONSTRAINT chemmapchemicals_pkey PRIMARY KEY (id),
    CONSTRAINT unique_dsstox_id UNIQUE (dsstox_id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chemicals
    OWNER to postgres;


GRANT USAGE, SELECT ON SEQUENCE chem_seq TO postgres;