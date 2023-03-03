-- Table: chem_toxexp_value

-- DROP TABLE chem_toxexp_value;

CREATE TABLE chem_toxexp_value
(
    dsstox_id text COLLATE pg_catalog."default",
    prop_value text[] COLLATE pg_catalog."default",
    CONSTRAINT dsstox_CATMOS UNIQUE (dsstox_id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chem_toxexp_value
    OWNER to postgres;