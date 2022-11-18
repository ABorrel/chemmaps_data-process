-- Table: public.chem_interference_prediction_name

-- DROP TABLE public.chem_interference_prediction_name;

CREATE TABLE public.chem_interference_prediction_name
(
    id numeric,
    name text COLLATE pg_catalog."default"
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chem_interference_prediction_name
    OWNER to chemmap_ntp;