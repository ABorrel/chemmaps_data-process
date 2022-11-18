-- Table: public.chem_toxexp_value

-- DROP TABLE public.chem_toxexp_value;

CREATE TABLE public.chem_toxexp_value
(
    dsstox_id text COLLATE pg_catalog."default",
    prop_value text[] COLLATE pg_catalog."default"
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chem_toxexp_value
    OWNER to postgres;