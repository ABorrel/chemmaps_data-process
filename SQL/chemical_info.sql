-- Table: public.chemical_info

-- DROP TABLE public.chemical_info;
create sequence chemmapchemical_info_id_seq;
CREATE TABLE public.chemical_info
(
    id bigint NOT NULL DEFAULT nextval('chemmapchemical_info_id_seq'::regclass),
    chem_id bigint NOT NULL,
    chemical_name text COLLATE pg_catalog."default",
    mol_weight numeric,
    source text COLLATE pg_catalog."default",
    source_id text COLLATE pg_catalog."default",
    epa_dashboard_id bigint,
    tox21 text COLLATE pg_catalog."default",
    pfas text COLLATE pg_catalog."default",
    sub_source text COLLATE pg_catalog."default",
    sub_source_id text COLLATE pg_catalog."default",
    CONSTRAINT chemmapchemical_info_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chemical_info
    OWNER to postgres;