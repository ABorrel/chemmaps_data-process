-- Table: public.ice_tox21

-- DROP TABLE ice_tox21;

CREATE TABLE ice_tox212
(
    dtxsid text COLLATE pg_catalog."default",
    data_sources text COLLATE pg_catalog."default",
    number_of_pubmed_articles numeric,
    pubchem_data_sources numeric,
    cpdat_count numeric,
    spid text COLLATE pg_catalog."default",
    qc_new_spid text COLLATE pg_catalog."default",
    chid text COLLATE pg_catalog."default",
    qc_grade text COLLATE pg_catalog."default",
    qc_grade_description text COLLATE pg_catalog."default",
    code text COLLATE pg_catalog."default",
    chnm text COLLATE pg_catalog."default",
    aeid text COLLATE pg_catalog."default",
    aenm text COLLATE pg_catalog."default",
    original_hitc numeric,
    new_hitc numeric,
    modl_tp numeric,
    modl_ga numeric,
    modl_acc numeric,
    modl text COLLATE pg_catalog."default",
    m4id text COLLATE pg_catalog."default",
    logc_min text COLLATE pg_catalog."default",
    logc_max text COLLATE pg_catalog."default",
    resp_unit text COLLATE pg_catalog."default",
    max_med numeric,
    max_med_conc numeric,
    coff numeric,
    qc_omit_src text COLLATE pg_catalog."default",
    new_hitc_flag text COLLATE pg_catalog."default",
    ac50 numeric,
    acc numeric,
    min_c_tested numeric,
    max_c_tested numeric
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE ice_tox212
    OWNER to postgres;