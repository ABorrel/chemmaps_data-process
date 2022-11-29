-- Table: chem_prop_drugbank_value

-- DROP TABLE chem_prop_drugbank_value;
CREATE TABLE chem_prop_drugbank_value
(
    drugbank_id text COLLATE pg_catalog."default",
    prop_value text[] COLLATE pg_catalog."default",
    CONSTRAINT unique_source_id UNIQUE (drugbank_id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chem_prop_drugbank_value
    OWNER to postgres;