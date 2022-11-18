-- Table: chem_prop_drugbank_name

-- DROP TABLE chem_prop_drugbank_name;
create sequence chemmap_drugbank_name_prop_id_seq;
CREATE TABLE chem_prop_drugbank_name
(
    id integer NOT NULL DEFAULT nextval('chemmap_drugbank_name_prop_id_seq'::regclass),
    name text COLLATE pg_catalog."default",
    CONSTRAINT chemmap_drugbank_name_prop_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chem_prop_drugbank_name
    OWNER to postgres;

GRANT USAGE, SELECT ON SEQUENCE chemmap_drugbank_name_prop_id_seq TO postgres;