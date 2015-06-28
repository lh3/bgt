## Querying Genotypes with Standard SQL

This document demonstrates that we can use three *conceptual* SQL tables for
variant annotations, sample phenotypes and genotypes, respectively, to achieve
flexible query of genotype data. The tables are *conceptual* in that they are
primarily used for constructing queries but not for storing real data as they
are not efficient enough for large-scale data.

We will use an example to show the structure of the three tables and how to
query data in standard SQL. You can find all the SQL statements in
[Appendix](#appendix) at the end of this document.

### An example

We define three tables as follows:
```sql
CREATE TABLE Variant (  -- variant annotations
  vid  TEXT, -- variant ID; unique
  gene TEXT  -- gene name; can be null
);
CREATE TABLE Sample (   -- sample phenotypes
  sid TEXT,  -- sample ID; unique
  age REAL   -- age of the sample
);
CREATE TABLE Genotype ( -- sample genotypes kept in a sparse way
  vid TEXT,  -- variant ID; joins Variant.vid
  sid TEXT,  -- sample ID;  joins Sample.sid
  gt  INT    -- number of non-reference alleles in the genotype
);

```
For demonstration, we assume the tables are filled with:
```
+-----+------+    +-----+-----+    +-----+-----+----+
| vid | gene |    | sid | age |    | vid | sid | gt |
+-----+------+    +-----+-----+    +-----+-----+----+
| V1  | TP53 |    | S1  | 31  |    | V1  | S1  | 2  |
| V2  | CDK2 |    | S2  | 27  |    | V1  | S2  | 1  |
+-----+------+    | S3  | 58  |    | V1  | S3  | 0  |
                  +-----+-----+    | V2  | S1  | 0  |
                                   | V2  | S2  | 1  |
                                   | V2  | S3  | 1  |
                                   +-----+-----+----+
```
With these three tables alone, we can do a variety of queries. For simple
examples without looking at genotypes:
```sql
-- get variants on gene TP53
SELECT vid          -- information to output
  FROM Variant      -- table queried
 WHERE gene="TP53"; -- condition
-- get samples younger than 40
SELECT sid
  FROM Sample
 WHERE age<40;
```
Due to the way `gt` is defined, sum of `gt` gives the count of non-reference allele.
We can use `GROUP BY` and [aggregate functions][aggfunc] to do computation on
genotypes.
```sql
-- get non-ref allele count for all variants
SELECT vid,SUM(gt) -- SUM of gt grouped by vid
  FROM Genotype
 GROUP BY vid;
-- for each variant, get the list of samples having the ALT
SELECT vid,GROUP_CONCAT(sid)
  FROM Genotype
 WHERE gt>0
 GROUP BY vid;
```
We can set conditions on aggregate functions using the `HAVING` clause:
```sql
-- get TP53 or CDC2 variants with allele count over 2
SELECT g.vid
  FROM Genotype g
 GROUP BY g.vid
 HAVING SUM(g.gt)>2;
```
We can apply table joining to string queries on multiple tables:
```sql
-- get TP53 variants that are common in young samples
SELECT g.vid
  FROM Variant v, Sample s, Genotype g
 WHERE v.gene="TP53" AND s.age<40
   AND v.vid=g.vid AND s.sid=g.sid -- table join
 GROUP BY g.vid
 HAVING SUM(g.gt)>2;
```
A table can have multiple instances and glued together. This gives more power:
```sql
-- get CDK2/TP53 variants that are common in young samples but rare in old
SELECT v.vid
  FROM Sample s1, Genotype g1,
       Sample s2, Genotype g2,
       Variant v
 WHERE s1.age<40 AND s1.sid=g1.sid
   AND s2.age>50 AND s2.sid=g2.sid
   AND (v.gene="CDK2" OR v.gene="TP53")
   AND v.vid=g1.vid AND v.vid=g2.vid
 GROUP BY v.vid
 HAVING SUM(g1.gt)>2 AND SUM(g2.gt)<1;
```

### Practical considerations

The `Genotype` table stores sample genotypes in a sparse way. Given tens of
thousands of whole-genome samples, this table will have trillions of rows.
A standard SQL engine is not efficient enough to handle these many rows and
to compute aggregate. We need to mimic the behavior of this `Genotype` table
with another data structure.

### Acknowledgement

The idea of SQL-based genotype query is inspired by [GEMINI][gemini] and
[GQT][gqt]. The conceptual SQL schema is influenced by [Mike Lin][mlin]'s SQL
example sent to the RefVar working group of GA4GH.

### <a name="appendix"></a>Appendix

To try the example, copy and paste the SQL statements in a file `test.sql` and execute them with
```sh
sqlite3 < test.sql
```
SQL statements:
```sql
DROP TABLE IF EXISTS Variant;
CREATE TABLE Variant (
	vid  TEXT,
	gene TEXT
);

DROP TABLE IF EXISTS Sample;
CREATE TABLE Sample (
	sid TEXT,
	age REAL
);

DROP TABLE IF EXISTS Genotype;
CREATE TABLE Genotype (
	vid TEXT,
	sid TEXT,
	gt  INT
);

INSERT INTO Variant VALUES ("V1", "TP53");
INSERT INTO Variant VALUES ("V2", "CDK2");

INSERT INTO Sample VALUES ("S1", 31);
INSERT INTO Sample VALUES ("S2", 27);
INSERT INTO Sample VALUES ("S3", 58);

INSERT INTO Genotype VALUES ("V1", "S1", 2);
INSERT INTO Genotype VALUES ("V1", "S2", 1);
INSERT INTO Genotype VALUES ("V1", "S3", 0);
INSERT INTO Genotype VALUES ("V2", "S1", 0);
INSERT INTO Genotype VALUES ("V2", "S2", 0);
INSERT INTO Genotype VALUES ("V2", "S3", 1);

SELECT vid FROM Variant WHERE gene="TP53";

SELECT sid FROM Sample WHERE age<40;

SELECT vid,SUM(gt)
  FROM Genotype
 GROUP BY vid;

SELECT vid,GROUP_CONCAT(sid)
  FROM Genotype
 WHERE gt>0
 GROUP BY vid;

SELECT vid
  FROM Genotype
 GROUP BY vid
 HAVING SUM(gt)>2;

SELECT g.vid
  FROM Variant v, Sample s, Genotype g
 WHERE v.gene="TP53" AND s.age<40
   AND v.vid=g.vid AND s.sid=g.sid -- table join
 GROUP BY g.vid
 HAVING SUM(g.gt)>2;

SELECT g1.vid,SUM(g1.gt),SUM(g2.gt)
  FROM Sample s1, Genotype g1,
       Sample s2, Genotype g2,
       Variant v
 WHERE s1.age<40 AND s1.sid=g1.sid
   AND s2.age>40 AND s2.sid=g2.sid
   AND (v.gene="CDK2" OR v.gene="TP53")
   AND v.vid=g1.vid AND v.vid=g2.vid
 GROUP BY v.vid
 HAVING SUM(g1.gt)>2 AND SUM(g2.gt)<1;

```

[gemini]: http://gemini.readthedocs.org/en/latest/
[gqt]: https://github.com/ryanlayer/gqt
[aggfunc]: https://en.wikipedia.org/wiki/Aggregate_function
[mlin]: http://www.mlin.net
