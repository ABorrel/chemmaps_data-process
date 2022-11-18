-- Table: public.bodymap_chemsum

-- DROP TABLE public.bodymap_chemsum;

CREATE SEQUENCE bodymap_genemap_1_id_seq;

CREATE TABLE public.bodymap_chemsum
(
    casn text COLLATE pg_catalog."default",
    ncct_qc_simple text COLLATE pg_catalog."default",
    smiles text COLLATE pg_catalog."default",
    name text COLLATE pg_catalog."default",
    dsstox_id text COLLATE pg_catalog."default",
    assays_tested text COLLATE pg_catalog."default",
    qc_t4_simple text COLLATE pg_catalog."default",
    qc_new_spid text COLLATE pg_catalog."default"
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.bodymap_chemsum
    OWNER to postgres;