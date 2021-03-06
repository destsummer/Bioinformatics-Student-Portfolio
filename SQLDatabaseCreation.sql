--This file includes the starting code to create a simple database and the input information for it.
--The database was created to offer beginning students in computer science a user friendly environment
--to learn about important figures in the field and their accomplishments.

--Final url for interactive use: http://pluto.hood.edu/~team6/project_home.php 

--Table input below:

Drop table if exists Computer_scientist;
Drop table if exists Schooling;
Drop table if exists Company;
Drop table if exists Field;
Drop table if exists Award;
Drop table if exists Developments;
Drop table if exists Applications;
Drop table if exists Spouse;
Drop table if exists Children;
Drop table if exists Attended;
Drop table if exists Founded;
Drop table if exists Works_in;
Drop table if exists Received;
Drop table if exists Discovered;
Drop table if exists Applied_to;

Create table Computer_Scientist(
	CS_ID float(8,2), 
	CS_name varchar(20) unique,
	DOB date not null,
	DOD date,
	CS_bio varchar(250),
	primary key (CS_ID));

Create table Schooling(
	University varchar(20),
	Country varchar(20),
	primary key(University));

Create table Company(
	Company_name varchar(20),
	Description varchar(250),
	primary key(Company_name));

Create table Field(
	Field_name varchar(20),
	Description varchar(250),
	primary key(Field_name));

Create table Award(
	Award_name varchar(20),
	primary key(Award_name));

Create table Developments(
	Development_name varchar(20),
	Description varchar(250),
	primary key(Development_name));

Create table Applications(
	Application_name varchar(20),
	Description varchar(250),
	primary key(Application_name));

Create table Spouse(
	CS_ID float(8,2),
	Spouse_ID float(8,2),
	Spouse _name varchar(20) UNIQUE,
	primary key(CS_ID,Spouse_ID),
	forgein key(CS_ID) references Computer _Scientist(CS_ID),
	ON DELETE CASCADE);

Create table Children(
	CS_ID float(8,2),
	Number_of_children float(8,2),
	primary key(CS_ID,Number_of_children),
	foreign key(CS_ID) references Computer_Scientist(CS_ID)
ON DELETE CASCADE);

Create table Founded(
	CS_ID float(8,2),
	Company_name varchar(20),
	CS_position varchar(20),
	primary key(CS_ID, Company_name),
	foreign key(CS_ID) references Computer_Scientist(CS_ID),
	foreign key(Company_name) references Company(Company_name));

Create table Attended(
	CS_ID float(8,2),
	University varchar(20),
	Year_of_graduation float(4),
	primary key(CS_ID, University),
	foreign key(CS_ID) references Computer_Scientist(CS_ID),
	foreign key(University) references Schooling(University));

Create table Received(
	CS_ID float(8,2),
	Award_name varchar(20),
	Year_received float(4),
	primary key(CS_ID, Award_name),
	foreign key(CS_ID) references Computer_Scientist(CS_ID),
	foreign key(Award_name) references Award(Award_name));	
	
Create table Works_in(
	CS_ID float(8,2),
	Field_name varchar (20), 
	primary key(CS_ID, Field_name),
	foreign key (CS_ID) references Computer_Scientist(CS_ID),
	foreign key (Field_name) references Field(Field_name));

Create table Discovered(
	CS_ID float(8,2),
	Development_name varchar(20),
	primary key(CS_ID, Development_name),
	foreign key(CS_ID) references Computer_Scientist(CS_ID),
	foreign key(Development_name) references Developments(Development_name));

Create table Applied_to(
	Development_name varchar(20),
	Application_name varchar(20),
	primary key(Development_name, Application_name),
	foreign key(Development_name) references Developments(Development_name),
	foreign key(Application_name) references Applications(Application_name));

--Insertion of data: This method was used to ensure that would be no confusion as to what was entered into each tuple.
--Simple one line code can be used to enter all data into each table.

insert into Computer_Scientist(CS_ID, CS_name, DOB, DOD, CS_bio) values ('0000001', 'Elon Musk', '1971-06-28', '/N’', 'Elon, form South Africa, started his career in his early 20s. Today, he is the founder of many companies that are known worldwide.' );
insert into Computer_Scientist(CS_ID, CS_name, DOB, DOD, CS_bio) values ('0000002', 'Brendan Eich', '1961-07-04', '/N', 'Brendan is known for Javascript which is used by HTML. He is also the cofounder of an internet browser called Mozilla.' );
insert into Computer_Scientist(CS_ID, CS_name, DOB, DOD, CS_bio) values ('0000003', 'Linus Torvalds', '1969-12-28', '/N', 'Linus is a Finnish mastermind behind the Linux operating system. He is also known for other important computer science implementations.');

insert into Schooling(University, Country) values ('University of Pennsylvania', 'US');
insert into Schooling(University, Country) values ('Stanford University', 'US');
insert into Schooling(University, Country) values ('Queens University', 'Canada');
insert into Schooling(University, Country) values ('Santa Clara University', 'US');
insert into Schooling(University, Country) values ('University of Illinois', 'US');
insert into Schooling(University, Country) values ('University of Helsinki', 'Finland');

insert into Company(Company_name, Description) values ('Tesla Motors', 'Company that formed in 2003 to provide affordable electric cars. Battery products and solar panels are also included.');
insert into Company(Company_name, Description) values ('Brave Software', 'Internet browser platform. Brave uses a Basic Attention Token which is a cryptocurrency used for raising money.');
insert into Company(Company_name, Description) values ('Paypal', 'Electronic commerce that allows for payments between parties through online transfers.');

insert into Field(Field_name, Description) values ('Programming Languages', 'Vocabulary and set of grammatical rules for instructing a computer or computing device to perform specific tasks. Typical are Java, C, C++.');
insert into Field(Field_name, Description) values ('Software Engineering', 'Develops information systems by designing, developing, and installing software solutions.');
insert into Field(Field_name, Description) values ('Systems', 'Building programs that use a lot of resources and profiling that resource usage. Systems work includes building operating systems, databases and more.');
insert into Field(Field_name, Description) values (‘Graphics’, ‘);
insert into Field(Field_name, Description) values (‘Networking’, ‘);
insert into Field(Field_name, Description) values (‘Hardware’, ‘);

insert into Award(Award_name) values ('Millennium Tech Prize');

insert into Developments(Development_name, Description) values ('Javascript','Originating at Netscape, this language combined Scheme, with and orientation of Self and syntax of Java. Javascript is used web browsers.');
insert into Developments(Development_name, Description) values ('Linux Kernel', 'Altered version of UNIX operating system, that is available free to users.');
insert into Developments(Development_name, Description) values ('Git', 'Version controlled system that allows for tracking files, source code management, and much more.');

insert into Applications(Application_name, Description) values ('Mozilla', 'Co-founded by Brendan Eich in 1998. This was created to allow implementation of source code from Netscape.');
insert into Applications(Application_name, Description) values ('Github', 'Online open source community for users to collaborate, share and discover.');
insert into Applications(Application_name, Description) values ('Linux', 'Open source operating system that uses Linux Kernel that directly links a systems hardware and resources.');

insert into Works_in(CS_ID, Field_name) values ('0000001', 'Software Engineering');
insert into Works_in(CS_ID, Field_name) values ('0000002', 'Programming Languages');
insert into Works_in(CS_ID, Field_name) values ('0000002', 'Software Engineering');
insert into Works_in(CS_ID, Field_name) values ('0000003', 'Systems');

insert into Discovered(CS_ID, Development_name) values ('0000002', 'Javascript');
insert into Discovered(CS_ID, Development_name) values ('0000003', 'Linux Kernel');
insert into Discovered(CS_ID, Development_name) values ('0000003', 'Git');

insert into Applied_to(Development_name, Application_name) values ('Javascript', 'Mozilla');
insert into Applied_to(Development_name, Application_name) values ('Linux Kernel', 'Linux');
insert into Applied_to(Development_name, Application_name) values ('Git', 'Github');

insert into Founded(CS_ID, Company_name, CS_position) values ('0000002', 'Brave Software', 'CEO');
insert into Founded(CS_ID, Company_name, CS_position) values ('0000001', 'Tesla Motors', 'CEO');
insert into Founded(CS_ID, Company_name, CS_position) values ('0000001', 'Paypal', 'CEO');

insert into Attended(CS_ID, University, Year_of_graduation) values ('0000002', 'Santa Clara University', '1983');
insert into Attended(CS_ID, University, Year_of_graduation) values ('0000002', 'University of Illinois', '1985');
insert into Attended(CS_ID, University, Year_of_graduation) values ('0000001', 'University of Pennsylvania', '1997');
insert into Attended(CS_ID, University, Year_of_graduation) values ('0000001', 'Stanford University', '1995');
insert into Attended(CS_ID, University, Year_of_graduation) values ('0000001', 'Queens University', '1992');
insert into Attended(CS_ID, University, Year_of_graduation) values ('0000003', 'University of Helsinki', '1988');

insert into Received(CS_ID, Award_name, Year_Received) values ('0000003', 'Millennium Tech Prize', '2012');

insert into Spouse(CS_ID, Spouse_ID, Spouse_name) values ('0000002', '0000021', 'Eleanor');
insert into Spouse(CS_ID, Spouse_ID, Spouse_name) values ('0000003', '0000031', 'Tove');
insert into Spouse(CS_ID, Spouse_ID, Spouse_name) values ('0000001', '0000011', 'Talulah');

insert into Children(CS_ID,Number_of_children) values ('0000001', '6');
insert into Children(CS_ID,Number_of_children) values ('0000002', '5');
insert into Children(CS_ID,Number_of_children) values ('0000003', '3');

