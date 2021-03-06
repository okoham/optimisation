/*
Die Ehepaare Richter, Engler, Passmann und Kaiser unterhalten sich darüber, in welchem
Land sie in diesem Jahr ihren Urlaub verbringen wollen. Wie lauten die Vor- und Nachnamen
der Eheleute, und in welchem Land wollen sie die schönste Zeit des Jahres verbringen?

Anton Kaiser und Regina Richter wollen nicht in Österreich Urlaub machen.

Frank und seine Frau verbringen ihren Urlaub ieutschland und Inge und ihr Mann fahren in die Schweiz.

Die Eheleute Kurt und Katja reisen in diesem Jahr nicht nach Spanien.

Inges Nachname lautet nicht Kaiser und Philipps Nachname nicht Engler.

Susanne freut sich schon sehr auf ihren Urlaub.
*/

include "alldifferent.mzn";

enum families = {richter, engler, passmann, kaiser};
enum girls = {inge, regina, katja, susanne};
enum boys = {anton, frank, kurt, philipp};
enum countries = {de, at, ch, es};

array[families] of var girls: wife;
array[families] of var boys: husband;
array[families] of var countries: destination;

constraint alldifferent(wife);
constraint alldifferent(husband);
constraint alldifferent(destination);

constraint husband[kaiser] == anton;
constraint wife[richter] == regina;
constraint destination[kaiser] != at;
constraint destination[richter] != at;

constraint forall (f in families) (husband[f] == frank -> destination[f] == de);
constraint forall (f in families) (wife[f] == inge -> destination[f] == ch);

constraint forall (f in families) (husband[f] == kurt <-> wife[f] = katja /\ destination[f] != es);

constraint wife[kaiser] != inge;
constraint husband[engler] != philipp;

solve satisfy;

output [ "family, wife, husband, destination \n" ];
output [ "\(f), \(wife[f]), \(husband[f]), \(destination[f])\n" | f in families ];