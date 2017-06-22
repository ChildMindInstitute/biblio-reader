import random, os, collections, math, manager as mg
checker_dir = mg.dir(os.path.join(mg.ROOT_PATH, 'checker_assigns'))
data = mg.get_data()
paragraphs = mg.get_paragraphs()


class Member(object):
    def __init__(self, name, path):
        """
        The Member class represents a reviewer of the articles
        :param name: Reviewer ID
        :param path: Path of reviewer text file
        """
        self.name = name
        self.path = '/'.join([path, name + '.txt'])
        if os.path.exists(self.path):
            file = open(self.path)
            if len(file.readlines()) < 4:
                self.articles = []
            else:
                for i, line in enumerate(file):
                    if i == 3:
                        self.articles = sorted([int(l) for l in str(line).strip().split(',')])
                        break
        else:
            self.articles = []
        self.written = list(self.articles)

    def __str__(self):
        """
        The reviwer's text file contains key information associated with their articles, including article Title,
        Publication type, authors, journal category, and key text snippets containing matched key words
        """
        return self.name + '\n\n\n' + ','.join([str(article) for article in sorted(self.articles)]) + '\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(key) + ': ' + str(data.loc[key, 'Title']) +
                        '\n' + str(data.loc[key, 'Authors']) + '\nPublication Type: ' +
                        data.loc[key, 'Journal Category'] + '\n\n' + '\n\n'.join(paragraph)
                        for key, paragraph in sorted(paragraphs.items()) if key in self.articles])

    def assign(self, assign_list, length):
        """
        Randomly assigns numbers from the list to the member's list of articles
        :param assign_list: The list of article reference numbers to assign to reviewers
        :param length: The size of numbers to assign
        """
        while length > 0 and not set(assign_list).issubset(set(self.articles)):
            assignment = random.choice(assign_list)
            if assignment not in self.articles:
                self.articles.append(assignment)
                assign.remove(assignment)
                length -= 1


class Assignment(object):
    def __init__(self, members=None, dir=checker_dir):
        """
        Creates members with member names and dir paths from the members list
        :param members: if None, looks through the directory for all text files contained and creates members based on
        those. Otherwise, creates new members
        :param dir: The path for each member
        """
        if members is None:
            self.members = [Member(member.replace('.txt', ''), dir) for member in os.listdir(dir) if
                            '.txt' in member]
        else:
            self.members = [Member(member, dir) for member in members]

    def __getitem__(self, item):
        if item in range(0, len(self.members)):
            return self.members[item]
        for member in self.members:
            if item == member.name:
                return member
        return None

    def assign(self, assignment, length=None):
        """
        Assigns each member an equal amount of articles from assignment
        :param assignment: The assignment of articles
        :param length: The amount of articles to assign
        """
        if length is None:
            length = len(assignment)
        length /= len(self.members)
        for member in self.members:
            member.assign(assignment, length)

    def write(self, new=None):
        """
        Creates text files for each member
        :param new: If not None, creates a separate new text file for each member if they already have one
        """
        for member in self.members:
            if new:
                member.articles = [article for article in member.articles if article not in member.written]
                with open(member.path.replace('.txt', '') + new + '.txt', 'w') as f:
                    f.write(str(member))
            else:
                with open(member.path, 'w') as f:
                    f.write(str(member))

    def double(self):
        """
        After assigning articles to members, reassigns articles to be double-checked (no member is assigned
        any article twice)
        """
        articles = []
        for member in self.members:
            articles += member.articles
        assignment = [item for item, count in collections.Counter(articles).items() if count == 1]
        member_length = math.ceil((len(articles) + len(assignment)) / len(self.members))
        revert = None
        while len(assignment) > 0:
            for member in self.members:
                if revert:
                     while revert > 0:
                         reverted = random.choice(member.articles)
                         if reverted not in assignment and reverted not in member.written:
                             assignment.append(reverted)
                             member.articles.remove(reverted)
                             revert -= 1
                     revert = None
                if len(member.articles) == member_length:
                    continue
                else:
                    member.assign(assignment, member_length - len(member.articles))
                    if len(member.articles) != member_length:
                        revert = len(assignment)

    def test(self, test_type=list([1, 2, 3])):
        """
        Tests to make sure no member has duplicate articles, and that the assignment has been assigned correctly
        :param test_type: Different tests for the assignment (default is all)
        """
        if 1 in test_type:
            articles = []
            for member in self.members:
                articles += member.articles
            print('Duplicate items:',
                  len([item for item, count in collections.Counter(articles).items() if count == 2]))
        if 2 in test_type:
            for member in self.members:
                print(member.name.capitalize() + "'s Length:", len(member.articles), '| Duplicates:',
                      len([item for item, count in collections.Counter(member.articles).items() if count > 1]))
        if 3 in test_type:
            for member in self.members:
                print(member.name.capitalize() + "'s New Articles:",
                      len([article for article in member.articles if article not in member.written]), '| Old:',
                      len(member.written))