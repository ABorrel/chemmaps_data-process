-- Table: public.bodymap_assay_mapping_new

-- DROP TABLE public.bodymap_assay_mapping_new;

CREATE TABLE public.bodymap_assay_mapping_new
(
    assay text COLLATE pg_catalog."default",
    gene text COLLATE pg_catalog."default",
    type_map text COLLATE pg_catalog."default",
    organ text COLLATE pg_catalog."default",
    system text COLLATE pg_catalog."default"
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.bodymap_assay_mapping_new
    OWNER to postgres;