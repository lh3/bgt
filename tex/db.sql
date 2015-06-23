CREATE TABLE Sample (
	sid TEXT,
	pop  TEXT
);

CREATE TABLE Genotype (
	vid TEXT,
	sid TEXT,
	gt  INT
);

INSERT INTO Sample VALUES ("s1", "CEU");
INSERT INTO Sample VALUES ("s2", "CEU");
INSERT INTO Sample VALUES ("s3", "AFR");

INSERT INTO Genotype VALUES ("v1", "s1", 2);
INSERT INTO Genotype VALUES ("v1", "s2", 1);
INSERT INTO Genotype VALUES ("v1", "s3", 0);
INSERT INTO Genotype VALUES ("v2", "s1", 0);
INSERT INTO Genotype VALUES ("v2", "s2", 0);
INSERT INTO Genotype VALUES ("v2", "s3", 1);

SELECT g1.vid,SUM(g1.gt),SUM(g2.gt)
  FROM Sample s1, Genotype g1,
       Sample s2, Genotype g2
 WHERE s1.pop="CEU" AND s1.sid=g1.sid
   AND s2.pop="AFR" AND s2.sid=g2.sid
   AND g1.vid=g2.vid
 GROUP BY g1.vid
 HAVING SUM(g1.gt)>2 AND SUM(g2.gt)<1;
