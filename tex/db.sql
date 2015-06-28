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

