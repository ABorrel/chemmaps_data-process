-- Table: public.chem_prop_drugbank_value

-- DROP TABLE public.chem_prop_drugbank_value;
create sequence chemmap_prop_value_id_seq;
CREATE TABLE public.chem_prop_drugbank_value
(
    id bigint NOT NULL DEFAULT nextval('chemmap_prop_value_id_seq'::regclass),
    source_id text COLLATE pg_catalog."default",
    prop_value text[] COLLATE pg_catalog."default",
    CONSTRAINT chemmap_prop_value_pkey PRIMARY KEY (id),
    CONSTRAINT unique_source_id UNIQUE (source_id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chem_prop_drugbank_value
    OWNER to postgres;