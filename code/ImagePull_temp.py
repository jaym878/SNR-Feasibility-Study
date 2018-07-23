class ImagePull:

    def arrays(self, pos):
        # I dont know whether I should pass "line" or load it here
        #line = np.genfromtxt(open('/home/murphyj/Desktop/Coding/SNR_list.csv', "r"), names=True, delimiter=',', dtype=None)
        self.list = []
        num = self.line.shape[0]
        for n in range(num):
            list.append(self.line[n][self.pos])
        self.list = np.asarray(self.list)
        return list

