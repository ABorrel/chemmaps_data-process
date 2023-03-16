-- Table: cHTS_assays

-- DROP TABLE cHTS_assays;

CREATE TABLE cHTS_assays
(
    assay text COLLATE pg_catalog."default" NOT NULL,
    assay_source text COLLATE pg_catalog."default",
    species text COLLATE pg_catalog."default",
    tissue text COLLATE pg_catalog."default",
    invitro_assay_format text COLLATE pg_catalog."default",
    gene text COLLATE pg_catalog."default",
    entrez_gene_id text COLLATE pg_catalog."default",
    mechanistic_target text COLLATE pg_catalog."default",
    CONSTRAINT assay_unique UNIQUE (assay)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE cHTS_assays
    OWNER to postgres;


-- Table cHTS_invitroDB
-- DROP TABLE cHTS_invitroDB;

CREATE TABLE chts_invitroDB_2
(
    recordid text COLLATE pg_catalog."default" NOT NULL,
    dsstox_id text COLLATE pg_catalog."default" NOT NULL,
    casn text COLLATE pg_catalog."default" NOT NULL,
    assay text COLLATE pg_catalog."default",
    curve_flag_description text COLLATE pg_catalog."default",
    chemical_qc_description  text COLLATE pg_catalog."default",
    test_range text COLLATE pg_catalog."default",
    test_range_unit text COLLATE pg_catalog."default",
    endpoint text COLLATE pg_catalog."default",
    response text COLLATE pg_catalog."default",
    response_unit text COLLATE pg_catalog."default",
    reference text COLLATE pg_catalog."default",
    
    CONSTRAINT recordid_unique_2 UNIQUE (recordid)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chts_invitroDB_2
    OWNER to postgres;



-- Table cHTS_chemicalQC
-- DROP TABLE cHTS_chemicalQC;

CREATE TABLE cHTS_chemicalQC
(
    ncgc_id text COLLATE pg_catalog."default" NOT NULL,
    tox21_id text COLLATE pg_catalog."default" NOT NULL,
    casn text COLLATE pg_catalog."default",
    dsstox_id text COLLATE pg_catalog."default",
    chem_type  text COLLATE pg_catalog."default",
    tox21_qc_t0 text COLLATE pg_catalog."default",
    tox21_qc_t4 text COLLATE pg_catalog."default",
    NICEATM_qc_t0 text COLLATE pg_catalog."default",
    NICEATM_qc_t4 text COLLATE pg_catalog."default",    
    NICEATM_qc_summary_call text COLLATE pg_catalog."default",

    CONSTRAINT ncgc_id UNIQUE (ncgc_id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE cHTS_chemicalQC
    OWNER to postgres;



-- remove nan in tables
UPDATE chts_assays set gene=null WHERE gene='nan';
UPDATE invitro_assay_format set gene=null WHERE invitro_assay_format='nan';

